#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of naivemapper.
#
#  naivemapper is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from collections import defaultdict
from CGRtools.RDFread import RDFread
from CGRtools.RDFwrite import RDFwrite
from CGRtools.CGRcore import CGRcore
from CGRtools.FEAR import FEAR
from core.fragger import Fragger
from core.bitstringen import Bitstringen
from core.wedding import Wedding
from sklearn.naive_bayes import BernoulliNB
from sklearn.utils.tests.test_linear_assignment import _hungarian
from random import shuffle

import argparse
import sys
import traceback
import networkx as nx
import pandas as pd
import numpy as np
import pickle


def mapper_core(**kwargs):
    c = 0
    fragger = Fragger(kwargs['min'], kwargs['max'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['b_length'])
    wedding = Wedding(kwargs['pairs'], kwargs['duplicate'])

    if kwargs['stage'] == 2:
        for _ in RDFread(open(kwargs['input'])).readdata():
            c += 1
            # подсчитаем кол-во реакций, чтоб в дальнейшем разбить их на части(4/5 обучение и 1/5 предсказание)

    if kwargs['stage'] == 1:
        # если стадия предсказания, то загружаем уже существующую модель в Наивный Бейсовский классификатор
        with open(kwargs['model'], 'rb') as train:
            nb = pickle.load(train)
    else:
        # если стадия обучения, то создаем новую модель
        nb = BernoulliNB(alpha=1.0, binarize=None)

    def worker(file):  # для души. увидим ошибки в RDF
        err = 0
        num = 0
        for num, data in enumerate(file.readdata(), start=1):
            if kwargs['debug'] or num % 10 == 1:
                print("reaction: %d" % num, file=sys.stderr)
            try:
                yield data
            except Exception:
                err += 1
                print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)

        print('%d from %d reactions processed' % (num - err, num), file=sys.stderr)

    def getXY(reaction):
        sub_graph = nx.union_all(reaction['substrats'])
        prod_graph = nx.union_all(reaction['products'])
        sub_frag = fragger.get(sub_graph)  # Словарь фрагментов от данного атома реагента
        prod_frag = fragger.get(prod_graph)  # Словарь фрагментов от данного атома продукта
        pairs, y_bit = wedding.get(sub_graph, prod_graph)  # Генерирует номера соотвествующих пар(а_реагента,а_продукта)
        x_bit = bitstring.get(sub_frag, prod_frag, pairs)  # Создают битовую строку для пары атомов реагента и продукта

        return x_bit, y_bit, pairs

    def mapping(pairs, y):
        tmp = defaultdict(dict)  # создается кв.матрица (кол-во_атомов_sub)Х(кол-во_атомов_prod)
        for (s_atom, p_atom), proba in zip(pairs, y):
            tmp[s_atom][p_atom] = - proba[1]  # если данная пара атомов не сгенерирована ранее в pairs то значение None
        prob_matrix = pd.DataFrame(tmp)
        indexes = _hungarian(prob_matrix.fillna(np.inf))
        p_reindex = prob_matrix.index.tolist()
        s_reindex = prob_matrix.columns.values.tolist()
        # Вычислить решение Манкрес проблемы классического назначения и возвращать индексы для спариваний
        # наименьшей стоимости. Возвращает 2D массив индексов
        _map = {p_reindex[x]: s_reindex[y] for x, y in indexes}  # mapping on products
        return _map

    def truth(f_test, f_pred):  # Проверка соответствия
        CGR = CGRcore(type='0', stereo=False, b_templates=None, balance=0, c_rules= None, e_rules=None)
        fear = FEAR()
        with open(f_pred) as predfile, open(f_test) as testfile:
            ok = 0
            nok = 0
            er = []
            for i, (pred, test) in enumerate(zip(RDFread(predfile).readdata(), RDFread(testfile).readdata()), start=1):
                predCGR = CGR.getCGR(pred)
                testCGR = CGR.getCGR(test)
                predHash = fear.getreactionhash(predCGR)
                testHash = fear.getreactionhash(testCGR)
                if predHash == testHash:
                    ok += 1
                else:
                    nok += 1
                    er.append(i)
            print("O'k: " + str(ok*100/(ok + nok)) + "%, \nNot O'k: " + str(nok*100/(ok + nok)) + "%")
            if kwargs['debug']:
                print(er)

    if kwargs['stage'] == 1:
        print("Testing set descriptor calculation")
        with open(kwargs['input']) as fr, open(kwargs['output'], 'w') as fw:  # Открываю входящий и исходящие файлы
            outputdata = RDFwrite(fw)
            for reaction in worker(RDFread(fr)):  # берем по 1 реакции из входящего файла
                x, _, pairs = getXY(reaction)  # генерирум битовую строку и пары соответствующих атомов
                y = nb.predict_log_proba(x)  # на основании модели из данной битовой строки получаем Y[True/False]
                _map = mapping(pairs, y)
                tmp = []
                for graph in reaction['products']:
                    tmp.append(nx.relabel_nodes(graph, _map, copy=True))
                reaction['products'] = tmp
                outputdata.writedata(reaction)
        truth(kwargs['input'], kwargs['output'])

    elif kwargs['stage'] == 2:
        indexes = list(range(c))  # генерации в список номеров(индексов) реакций
        folds = kwargs['folds']
        repeat = kwargs['repeat']

        for r in range(repeat):
            print('Repeat ', r+1, '/', repeat)
            shuffle(indexes)  # перемешиваем индексы реакций

            print("Training set descriptor calculation")
            with open('cross_v/mapping'+str(r)+'.rdf', 'w') as fw:
                test_file = RDFwrite(fw)
                for x in range(folds):  # генерация фолдов(блоков) кросс-валидации
                    print('Fold ', x+1, '/', folds)
                    test = indexes[x::folds]  # для предсказания выбираются каждая пятая реакция (из перемешенного типа)
                    with open(kwargs['input']) as fr:
                        for num, reaction in enumerate(worker(RDFread(fr))):
                            if num in test:  # если номер реакции во входящем файле совпал с номером тестового набора
                                test_file.writedata(reaction)  # мы записываем реакцию в файл для маппинга
                            else:  # если номер реакции во входящем файле НЕ совпал с номером тестового набора
                                x, y, pairs = getXY(reaction)  # генерируем битовые строки Х и строку У
                                nb.partial_fit(x, y, classes=pd.Series([False, True]))  #Обучаем нашу модель на основании x-bit и y-bit

            print("Testing set descriptor calculation")
            with open('cross_v/mapping'+str(r)+'.rdf') as fr, open('cross_v/output'+str(r)+'.rdf', 'w') as fw:
                outputdata = RDFwrite(fw)
                for reaction in worker(RDFread(fr)):
                    x, _, pairs = getXY(reaction)  # генерирум битовую строку
                    y = nb.predict_log_proba(x)  # на основании модели из данной битовой строки получаем строку Y[T/F]
                    _map = mapping(pairs, y)
                    tmp = []
                    for graph in reaction['products']:
                        tmp.append(nx.relabel_nodes(graph, _map, copy=True))
                        # на основании обученной модели перемаппливаются атомы продукта
                    reaction['products'] = tmp
                    outputdata.writedata(reaction)
            truth('cross_v/mapping'+str(r)+'.rdf', 'cross_v/output'+str(r)+'.rdf')  # проверка предсказанных данных

    else:
        print("Training set descriptor calculation")
        with open(kwargs['input']) as fr:
            for reaction in worker(RDFread(fr)):
                x, y, _ = getXY(reaction)  # для каждой реакции генерируем битовые строки Х и строку У
                nb.partial_fit(x, y, classes=pd.Series([False, True]))  # обучаем нашу модель на основании x-bit и y-bit

        with open(kwargs['model'], 'wb') as train:
            pickle.dump(nb, train)  # записываем нашу обученную модель в файл


def parse_args():
    parser = argparse.ArgumentParser(description="", epilog="", prog='naivemapper')
    parser.add_argument("--version", "-v", action="version", version="0.1", default=False)
    parser.add_argument("--input", "-i", default="input.rdf", type=str, help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=str, help="RDF outputfile")
    parser.add_argument("--model", "-n", default="model.dat", type=str, help="Model file")
    parser.add_argument("--min", "-m", type=int, default=1, help="minimal fragments length")
    parser.add_argument("--max", "-M", type=int, default=8, help="maximal fragments length")
    parser.add_argument("--bitstring", "-b", type=int, default=2,
                        help="type of united bitstring for atomic bitstrings A and B: 0-A*B, 1-A+B+A*B, 2-A!B+B!A+A*B")
    parser.add_argument("--b_length", "-l", type=int, default=2048, help="lenght of bitstrings")
    parser.add_argument("--pairs", "-p", type=int, default=0,
                        help="type of united respective pairs: 0-dumb, 1-equivalent")
    parser.add_argument("--duplicate", "-d", type=int, default=1,
                        help="type of united respective pairs: 0-does duplicate, 1-has a duplicate")
    parser.add_argument("--stage", "-s", type=int, default=2,
                        help="type of stage of the process: 0 - learning, 1 - prediction, 2 - cross_validation")
    parser.add_argument("--folds", "-N", type=int, default=5,
                        help="the number of partitions of the data: sets to training and test")
    parser.add_argument("--repeat", "-r", type=int, default=2,
                        help="the number of repetitions of the cross-validation procedure")
    parser.add_argument("--debug", action='store_true', help="debug mod")
    return parser


def main():
    parser = parse_args()
    args = parser.parse_args()
    mapper_core(**vars(args))


if __name__ == '__main__':
    main()
