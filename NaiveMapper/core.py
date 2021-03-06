import sys
import traceback
from collections import defaultdict
from timeit import default_timer as timer

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.utils.tests.test_linear_assignment import _hungarian

from CGRtools.files.RDFrw import RDFread
from CGRtools.preparer import CGRpreparer
from .ValFix import ValFix
from .SymFix import symFix


def worker(file, debug=False):  # для души. увидим ошибки в RDF
    err, num = 0, 0
    for num, data in enumerate(file, start=1):
        if debug or num % 10 == 1:
            # При условии debug выводит информацию о кажой реакции, иначе о каждой 10-ой
            print("reaction: %d" % num, file=sys.stderr)
        try:
            yield data
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)

    print('%d from %d reactions processed' % (num - err, num), file=sys.stderr)


def getXY(reaction, fragger, pairwise, bitstring, chunk=None):
    # Общий граф молекул реагента и молекул продукта, соответственно
    sub_graph, prod_graph = nx.union_all(reaction['substrats']), nx.union_all(reaction['products'])
    # Словарь фрагментов, сгенерированных от каждого атома в графе реагента и продукта, соответственно
    sub_frag, prod_frag = fragger.get(sub_graph), fragger.get(prod_graph)
    # Генерирует список пар(атом_реагента,атом_продукта) и соответствующий ему список значений ААО (True/False)
    pairs, y_bit = pairwise.get(sub_graph, prod_graph)
    # x-bit - битовую строку (хешированный молекулярный отпечаток) для каждой пары атомов
    if chunk:
        # Необходимое разбиение процесса, из-за массивных битовых строк (kwargs['length']>100000)
        for i in range(0, len(pairs), chunk):
            yield bitstring.get(sub_frag, prod_frag, pairs[i:i + chunk]), y_bit[i:i + chunk], pairs[i:i + chunk]
    else:
        x_bit = bitstring.get(sub_frag, prod_frag, pairs)
        yield x_bit, y_bit, pairs


def mapping(pairs, y, prod_graph, sub_graph, weights=[0.5, 0.5, 1.5]):
    valFix = ValFix()
    tmp = defaultdict(dict)  # создается кв.матрица (кол-во_атомов_реагента)Х(кол-во_атомов_продукта)
    for (s_atom, p_atom), proba in zip(pairs, y):
        tmp[s_atom][p_atom] = - proba[1]  # если данная пара атомов не сгенерирована ранее в pairs то значение None
    prob_matrix = pd.DataFrame(tmp).fillna(np.inf)  # заменяем значение None на +бесконечность (np.inf)

    # Вычисление решения Манкрес, который возвращает 2D массив - индексы для спариваний с наименьшей стоимостью.
    indexes = _hungarian(prob_matrix)

    p_reindex = prob_matrix.index.tolist()  # наименование строк, отвечающие за нумерацию атомов продукта
    s_reindex = prob_matrix.columns.values.tolist()  # наименование столбцов,отвечающие за нумерацию атомов реагента

    _m = {p_reindex[p]: s_reindex[s] for p, s in indexes}  # словарь со значениями атомного отображения
    _map = symFix(sub_graph, prod_graph, _m, prob_matrix)
    # пересмотр решения Манкреса (поиск в глубину по графу продукта)
    _map2 = valFix.getMap(sub_graph, prod_graph, _map, prob_matrix, weights)
    return _map2


def truth(f_test, f_pred, ok, nok, er, debug=False):  # Проверка соответствия
    cgr = CGRpreparer()
    with open(f_pred, encoding='cp1251') as predfile, open(f_test, encoding='cp1251') as testfile:
        for i, (pred, test) in enumerate(zip(RDFread(predfile).read(), RDFread(testfile).read()), start=1):
            predHash = cgr.getCGR(pred).get_signature_hash()
            testHash = cgr.getCGR(test).get_signature_hash()
            if predHash == testHash:
                ok += 1
            else:
                nok += 1
                er.append(i)

        print("Percentage\n\tO'k: %0.5f , \nNot O'k: %0.5f" % ((ok*100/(ok + nok)), (nok*100/(ok + nok))))
        if debug:
            print("The number of errors {} in reactions\n{}".format(len(er), er))

    return ok, nok
