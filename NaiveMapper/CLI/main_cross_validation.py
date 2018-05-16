#!/usr/bin/env python3
#
from random import shuffle
from sklearn.naive_bayes import BernoulliNB
from sklearn.neural_network import MLPClassifier
from ..keras_mlp_class.keras_mlp import Keras_MLP
import networkx as nx
import numpy as np
import pandas as pd

from CGRtools.files.RDFrw import RDFread, RDFwrite
from ..bitstringen import Bitstringen
from ..core import getXY, mapping, truth, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def remap(graphs, maps):
    tmp = []
    for graph in graphs:
        tmp.append(graph.remap({k: v for k, v in maps.items() if k in graph}, copy=True))
    return tmp


def cross_validation_core(**kwargs):
    fragger = Fragger(f_type=kwargs['fragment_type'], _min=kwargs['min'], _max=kwargs['max'], _deep=kwargs['deep'],
                      _min2=kwargs['min2'], _max2=kwargs['max2'], _fuzzy=kwargs['fuzzy'])
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
            file_1 = '{}/mapping{}{}.rdf'.format(kwargs['output'], r, n)

            if kwargs['type_model'] is 'nb':
                m = BernoulliNB(alpha=1.0, binarize=None)  # Создаем новую модель Наивного Бейсовского классификатора
            elif kwargs['type_model'] is 'mlp':
                hls, a, s, alpha = tuple(kwargs['mlp_hls']), kwargs['mlp_a'], kwargs['mlp_s'], kwargs['mlp_alpha']
                bs, lr, mi, es = kwargs['mlp_bs'], kwargs['mlp_lr'], kwargs['mlp_mi'], kwargs['mlp_es']
                m = MLPClassifier(hidden_layer_sizes=hls, activation=a, solver=s, alpha=alpha,
                                  batch_size=bs, learning_rate=lr, max_iter=mi, early_stopping=es)
            else:
                m = Keras_MLP(task="classification_2", layer_sizes=(400, 400), activations=['relu', 'relu'],
                              dropout="Auto", alpha=0.00001 * (2 ** 1), batch_size=200, learning_rate_init=0.001,
                              epochs=1,
                              shuffle=True, loss_function="mse", metrics=['mse'], verbose=1,
                              early_stopping=False, optimizer_name="adam", lr=0.001, beta_1=0.9, beta_2=0.999,
                              epsilon=1e-08)

            with open(file_1, 'w') as fw:  # with open(kwargs['input']) as fr, open(file_1, 'w') as fw:
                test_file = RDFwrite(fw)
                test = indexes[n::folds]  # для предсказания выбираются каждая N-ая реакция из списка
                train = [i for i in indexes if i not in test]

                X_train, Y_train = [], []
                for num, reaction in enumerate(worker(RDFread(open(kwargs['input'])), False)):
                    if num in test:  # Если номер рассматриваемой реакции совпал с номером тестового набора, то ...
                        # записываем её в файл для предсказания.
                        test_file.write(reaction)
                    else:  # если номер рассматриваем реакции НЕ совпал с номером тестового набора, то ...
                        for x, y, _ in getXY(reaction, fragger, pairwise1, bitstring, kwargs['chunk']):
                            X_train.append(x)
                            Y_train.append(y)
                            if num+1 % kwargs['batch_chunk'] == 0 or num == max(train):
                                """Обучаем нашу модель на основании: 
                                    - битовыx строк дескрипторов(Х), 
                                    - строк значений ААО(Y)"""
                                m.partial_fit(pd.concat(X_train, ignore_index=True),
                                              pd.concat(Y_train, ignore_index=True),
                                              classes=pd.Series([False, True]))
                                X_train.clear(), Y_train.clear()

            print("Testing set descriptor calculation")
            file_2 = '{}/output{}{}.rdf'.format(kwargs['output'], r, n)  # Контрольная выборка, с предсказанными ААО
            with open(file_1) as fr, open(file_2, 'w') as fw:
                output = RDFwrite(fw)
                for reaction in worker(RDFread(fr), kwargs['debug']):  # берем по 1 реакции из файла тестовго набора
                    p_graph = nx.union_all(reaction['products'])
                    p_nodes = [k + max(p_graph.nodes()) for k in p_graph.nodes()]
                    shuffle(p_nodes)
                    reaction.products._MindfulList__data = remap(reaction['products'],
                                                                 {k: v for k, v in zip(p_graph.nodes(), p_nodes)})

                    y, pairs = [], []
                    for x, _, drop_pairs in getXY(reaction, fragger, pairwise2, bitstring, kwargs['chunk']):
                        '''
                        на основании сгенерированного набора битовых строк из обученой модели выводятся значения 
                        логорифмов вероятностей проецирования (отображения) атомов
                        '''
                        pairs.extend(drop_pairs)
                        if kwargs['type_model'] is 'keras':
                            yi = m.predict(x.values)
                            y.extend([(0, p2[0]) for p2 in np.log(np.where(yi > 0, yi, 1e-20))])
                        else:
                            y.extend([yy for yy in m.predict_log_proba(x)])

                    _map = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
                    # на основании обученной модели перемаппливаются атомы продукта
                    reaction.products._MindfulList__data = remap(reaction['products'], _map)
                    output.write(reaction)  # запись реакции в исходящий файл

            ok, nok = truth(file_1, file_2, ok, nok, errors[n], kwargs['debug'])  # проверка предсказанных данных
