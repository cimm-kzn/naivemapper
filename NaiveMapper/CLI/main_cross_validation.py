import pickle

import networkx as nx
import pandas as pd
from sklearn.naive_bayes import BernoulliNB

from CGRtools.files.RDFrw import RDFread, RDFwrite
from ..bitstringen import Bitstringen
from ..core import getXY, mapping, truth, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def cross_validation_core(**kwargs):
    fragger = Fragger(kwargs['min'], kwargs['max'], kwargs['deep'], kwargs['fragment_type'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise1, pairwise2 = Pairwise(kwargs['pairs'], kwargs['duplicate']), Pairwise(0, 0)
    # pairwise2 - при предсказании, алгоритм Моргана(для выделения групп сим./экв. атомов) не применяется
    c = 0
    chunk = 200 if kwargs['bitstring'] is 'kron' else None

    # подсчитаем кол-во реакций, чтоб в дальнейшем разбить их на части(N-1/N обучение и 1/N предсказание)
    for _ in RDFread(open(kwargs['input'])).read():
        c += 1

    indexes = list(range(c))  # генерации в список номеров(индексов) реакций
    # число разбиений на блоки тестового набора и кол-во повторений процедуры валидации
    folds, repeat = kwargs['folds'], kwargs['repeat']

    for r in range(repeat):  # Генерация повторения процедуры валидации
        print('Repeat ', r + 1, '/', repeat)
        ok, nok = 0, 0
        errors = [[] for _ in range(folds)]
        # shuffle(indexes)  # перемешиваем индексы реакций

        for n in range(folds):  # генерация блоков(фолдов/разбиений) кросс-валидации
            print('Fold ', n + 1, '/', folds)
            print("Training set descriptor calculation")
            # Контрольная выборка, для оценки предсказательной способности
            file_1 = 'cross_v/mapping' + str(r) + str(n) + '.rdf'
            nb = BernoulliNB(alpha=1.0, binarize=None)  # Создаем новую модель Наивного Бейсовского классификатора

            with open(kwargs['input']) as fr, open(file_1, 'w') as fw:
                test_file = RDFwrite(fw)
                test = indexes[n::folds]  # для предсказания выбираются каждая N-ая реакция(из перемешенного списка)

                for num, reaction in enumerate(worker(RDFread(fr))):
                    if num in test:  # если номер рассматриваемой реакции совпал с номером тестового набора, то ...
                        # записываем её в файл для предсказания
                        test_file.write(reaction)
                    else:  # если номер рассматриваем реакции НЕ совпал с номером тестового набора, то ...
                        for x, y, _ in getXY(reaction, fragger, pairwise1, bitstring, chunk):
                            # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)
                            nb.partial_fit(x, y, classes=pd.Series([False, True]))

            print("Testing set descriptor calculation")
            file_2 = 'cross_v/output' + str(r) + str(n) + '.rdf'  # Контрольная выборка, с предсказанными ААО
            with open(file_1) as fr, open(file_2, 'w') as fw:
                output = RDFwrite(fw)
                for reaction in worker(RDFread(fr)):  # берем по 1 реакции из файла тестовго набора
                    y, pairs = [], []
                    for x, _, drop_pairs in getXY(reaction, fragger, pairwise2, bitstring, chunk):
                        '''на основании сгенерированного набора битовых строк дескрипторов
                        из обученой модели выводятся значения лагорифмов вероятностей проецирования (отображения) атомов'''
                        pairs.extend(drop_pairs)
                        y.extend([yy for yy in nb.predict_log_proba(x)])

                    _map = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
                    tmp = []
                    for graph in reaction['products']:
                        tmp.append(nx.relabel_nodes(graph, _map, copy=True))
                        # на основании обученной модели перемаппливаются атомы продукта
                    reaction['products'] = tmp
                    output.write(reaction)  # запись реакции в исходящий файл
            ok, nok = truth(file_1, file_2, ok, nok, errors[n])  # проверка предсказанных данных
