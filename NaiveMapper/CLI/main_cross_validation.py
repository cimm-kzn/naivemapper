from random import shuffle
from sklearn.naive_bayes import BernoulliNB
import networkx as nx
import pandas as pd

from CGRtools.files.RDFrw import RDFread, RDFwrite
from ..bitstringen import Bitstringen
from ..core import getXY, mapping, truth, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def remap(graphs, maps):
    tmp = []
    for graph in graphs:
        tmp.append(graph.remap(maps, copy=True))
    return tmp


def cross_validation_core(**kwargs):
    fragger = Fragger(kwargs['min'], kwargs['max'], kwargs['deep'], kwargs['fragment_type'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise1, pairwise2 = Pairwise(kwargs['pairs'], kwargs['duplicate']), Pairwise(0, False)
    # pairwise2 - при предсказании, алгоритм Моргана(для выделения групп сим./экв. атомов) никогда не применяется

    # подсчитаем кол-во реакций, чтоб в дальнейшем разбить их на части(N-1/N обучение и 1/N предсказание)
    c = 0
    for _ in RDFread(open(kwargs['input'])).read():
        c += 1

    # генерации в список номеров(индексов) реакций
    indexes = list(range(c))
    # число разбиений на блоки тестового набора и кол-во повторений процедуры валидации
    folds, repeat = kwargs['fold'], kwargs['repeat']

    for r in range(repeat):  # Генерация повторения процедуры валидации
        print('Repeat ', r + 1, '/', repeat)
        ok, nok = 0, 0
        errors = [[] for _ in range(folds)]
        # shuffle(indexes)  # перемешиваем индексы реакций

        for n in range(folds):  # генерация блоков (фолдов/разбиений) кросс-валидации
            print('Fold ', n + 1, '/', folds)
            print("Training set descriptor calculation")
            # Контрольная выборка, для оценки предсказательной способности
            file_1 = 'cross_v/mapping{}{}.rdf'.format(r, n)
            nb = BernoulliNB(alpha=1.0, binarize=None)  # Создаем новую модель Наивного Бейсовского классификатора

            with open(file_1, 'w') as fw:  # with open(kwargs['input']) as fr, open(file_1, 'w') as fw:
                test_file = RDFwrite(fw)
                test = indexes[n::folds]  # для предсказания выбираются каждая N-ая реакция(из перемешенного списка)

                for num, reaction in enumerate(worker(RDFread(open(kwargs['input'])), kwargs['debug'])):
                    if num in test:  # Если номер рассматриваемой реакции совпал с номером тестового набора, то ...
                        # записываем её в файл для предсказания.
                        test_file.write(reaction)
                    else:  # если номер рассматриваем реакции НЕ совпал с номером тестового набора, то ...
                        for x, y, _ in getXY(reaction, fragger, pairwise1, bitstring, kwargs['chunk']):
                            # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)
                            nb.partial_fit(x, y, classes=pd.Series([False, True]))

            print("Testing set descriptor calculation")
            file_2 = 'cross_v/output{}{}.rdf'.format(r, n)  # Контрольная выборка, с предсказанными ААО
            with open(file_1) as fr, open(file_2, 'w') as fw:
                output = RDFwrite(fw)
                for reaction in worker(RDFread(fr), kwargs['debug']):  # берем по 1 реакции из файла тестовго набора
                    p_graph = nx.union_all(reaction['products'])
                    reaction['products'] = remap(reaction['products'],
                                                 {k: k + max(p_graph.nodes()) for k in p_graph.nodes()})
                    y, pairs = [], []
                    for x, _, drop_pairs in getXY(reaction, fragger, pairwise2, bitstring, kwargs['chunk']):
                        '''
                        на основании сгенерированного набора битовых строк из обученой модели выводятся значения 
                        логорифмов вероятностей проецирования (отображения) атомов
                        '''
                        pairs.extend(drop_pairs)
                        y.extend([yy for yy in nb.predict_log_proba(x)])

                    _map = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
                    # на основании обученной модели перемаппливаются атомы продукта
                    reaction['products'] = remap(reaction['products'], _map)

                    output.write(reaction)  # запись реакции в исходящий файл
            ok, nok = truth(file_1, file_2, ok, nok, errors[n], kwargs['debug'])  # проверка предсказанных данных
