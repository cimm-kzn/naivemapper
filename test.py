# from CGRtools.strings import hash_cgr_string, get_morgan, get_cgr_string
# from CGRtools.containers import MoleculeContainer
# from CGRtools.preparer import CGRcombo
from CGRtools.files.RDFrw import RDFread, RDFwrite
from CGRtools.preparer import CGRpreparer
from NaiveMapper.SymFix import get_map_dfs
from NaiveMapper.ValFix import DFSdb
from NaiveMapper.bitstringen import Bitstringen
from NaiveMapper.core import getXY, truth, worker
from NaiveMapper.fragger import Fragger
from NaiveMapper.pairwise import Pairwise
from collections import defaultdict
from pprint import pprint
from sklearn.utils.tests.test_linear_assignment import _hungarian
from timeit import default_timer as timer

import hashlib as hs
import sys
import pickle
import pandas as pd
import networkx as nx
import numpy as np


# Сравнение ААО для истинных и предсказанных знвчений (ч.з сравнение хешей реакций)
""""""
cgr = CGRpreparer()
with open(sys.argv[1], encoding='cp1251') as fr, open(sys.argv[2], encoding='cp1251') as fw:
    ok, nok = 0, 0
    er = []

    for i, (pred, test) in enumerate(zip(RDFread(fw).read(), RDFread(fr).read()), start=1):
        predHash = cgr.getCGR(pred).get_signature_hash()
        testHash = cgr.getCGR(test).get_signature_hash()
        if predHash == testHash:
            ok += 1
        else:
            nok += 1
            er.append(i)

    print("Percentage\n\tO'k: %0.5f , \nNot O'k: %0.5f" % ((ok*100/(ok + nok)), (nok*100/(ok + nok))))
    print(len(er), '\n', er)


# Выявляем реакционные центры в ошибочных реакциях, и классифицируем их
"""
l = [23, 25, 37, 48, 58, 86]
cgr = CGRpreparer()
with open(sys.argv[1], encoding='cp1251') as fr:
    dictHash = defaultdict(list)
    for i, reac in enumerate(RDFread(fr), start=1):

        if i in l:
            '''
            rCGR = cgr.getCGR(reac)
            cgrRCenter = rCGR.get_environment(rCGR.get_center_atoms(), dante=False, deep=1)
            strRCenter = cgrRCenter.get_fear(cgrRCenter.get_morgan())
            hsRCenter = int(hs.md5(strRCenter.encode()).hexdigest(), 16)
            # hsRCenter = hash_cgr_string(strRCenter)
            dictHash[hsRCenter].append(i)
            '''
            rCGR = cgr.getCGR(reac)
            cgrRC = rCGR.get_environment(rCGR.get_center_atoms(), dante=False, deep=1)

            if len(set(cgr.split(cgrRC))) > 1:
                dictHash['more_1_RC'].append(i)
            else:
                dictHash['1_RC'].append(i)
            
            for rc in set(cgr.split(cgrRC)):
                rcStr = rc.get_fear(rc.get_morgan())
                dictHash[rcStr].append(i)
            '''
            strRC = cgrRC.get_fear(cgrRC.get_morgan())
            dictHash[strRC].append(i)
            '''
for k, v in sorted(dictHash.items(), key=lambda x: (len(x[1]), x[1][0]), reverse=False):
    print('"{}": {},'.format(k, v))
"""

# Выявление уникальных реакций
"""
with open(sys.argv[1], encoding='cp1251') as fr, open(sys.argv[2], "w") as fw:
    uniqueHash = {}
    print('Seeking unique items')
    for num, reaction in enumerate(RDFread(fr).read(), start=1):
        rHash = reaction.get_fear_hash()
        uniqueHash[rHash] = reaction

    print('Record file')
    outputdata = RDFwrite(fw)
    for v in uniqueHash.values():
        outputdata.write(v)

    print(len(uniqueHash), ' unique reactions')
"""

# Предсказание (проверка работы алгоритма DFS2)
"""
def remap(graphs, maps):
    tmp = []
    for graph in graphs:
        tmp.append(graph.remap(maps, copy=True))
    return tmp


dfs2 = DFSdb()
model = pickle.load(open(sys.argv[2], 'rb'))
fragger = Fragger(0, 3, 8, 2)
bitstring = Bitstringen(0, 2048, False)
pairwise = Pairwise(0, False)

num = 0
total_time = dict.fromkeys(["Frag+Bitstr", "Munkres", "Predict", "Remap1", "Remap2", "DFS1", "DFS2", "All"], 0)
time_ = timer()
with open(sys.argv[1], encoding='cp1251') as fr, open(sys.argv[3], 'w') as fw1:  # , open(sys.argv[4], 'w') as fw2:
    # out1, out2 = RDFwrite(fw1), RDFwrite(fw2)  # out = RDFwrite(fw1)
    out = RDFwrite(fw1)
    for num, reaction in enumerate(worker(RDFread(fr), True), start=1):
        y, pairs = [], []
        start_1 = timer()
        time_pr = 0
        for x, _, drop_pairs in getXY(reaction, fragger, pairwise, bitstring, False):
            pairs.extend(drop_pairs)
            start_2 = timer()
            y.extend(model.predict_log_proba(x))
            end_2 = timer()
            total_time["Predict"] += (end_2 - start_2)
            time_pr += (end_2 - start_2)
        total_time["Frag+Bitstr"] += (timer() - start_1 - time_pr)
        tmp = defaultdict(dict)
        for (s, p), proba in zip(pairs, y):
            tmp[s][p] = -proba[1]
        matrix = pd.DataFrame(tmp).fillna(np.inf)

        p_in, s_in = matrix.index.tolist(), matrix.columns.values.tolist()
        subG, prodG = nx.union_all(reaction['substrats']), nx.union_all(reaction['products'])

        start_ = timer()
        ind = _hungarian(matrix)
        total_time["Munkres"] += (timer() - start_)
        _m = {p_in[p]: s_in[s] for p, s in ind}  # словарь со значениями атомного отображения
        # print("Munckris:\t{}".format('_'.join(str(i) for _, i in sorted({s: p for p, s in _m.items()}.items()))))

        start_ = timer()
        _map = get_map_dfs(subG, prodG, _m)
        total_time["DFS1"] += (timer() - start_)
        '''
        start_ = timer()
        reaction['products'] = remap(reaction['products'], _map)
        total_time["Remap1"] += (timer() - start_)
        out.write(reaction)
        # пересмотр решения Манкреса (поиск в глубину по графу продукта)
        '''
        # if num != 13:
        start_ = timer()
        _map2 = dfs2.getMap(subG, prodG, _map, matrix)  # _m, matrix)
        total_time["DFS2"] += (timer() - start_)
        start_ = timer()
        reaction['products'] = remap(reaction['products'], _map2)
        total_time["Remap2"] += (timer() - start_)
        out.write(reaction)  # out2.write(reaction)

total_time["All"] += timer() - time_
start_ = timer()
_, _ = truth(sys.argv[1], sys.argv[3], 0, 0, [], True)
print("Проверка для {} реакций длилось {}".format(num, timer()-start_))
# _, _ = truth(sys.argv[1], sys.argv[4], 0, 0, [], True)  # Проверка соответствия
pprint(total_time)
"""

# Разбиваем rdf-файл на набор rxn-файлов
"""
file, fileOut = sys.argv[1], sys.argv[2]  # 'test/Ester/ester-zero.rdf'
with open(file, encoding='cp1251') as f1:
    tf, c, strLines = False, 0, str()
    for line in f1:
        if line.startswith('$RXN'):
            c += 1
            print(c, ' is done!')
            tf = True
        elif line.startswith('$DTYPE'):
            tf = False
        elif line.startswith('$RFMT') and c:
            out = "{}{}.rxn".format(fileOut, c)
            with open(out, 'w') as f2:
                f2.write(strLines)
            strLines = str()
        if tf:
            strLines += line
"""

# Выявление реакционных центров, запись их в отдельный файл.
# А так же создание rdf-файла с реакциями, имеющих общий реакционный центр.
"""
cgr = CGRpreparer()
reactions = RDFread(open(sys.argv[1])).read()
dictUnq = defaultdict(list)

with open(sys.argv[2], 'w') as fw:
    out = RDFwrite(fw)
    dictUnq2 = defaultdict(list)
    for num, reac in enumerate(reactions, start=1):
        if num%10 == 1: 
            print('Reactions %d is done' % (num))
        rCGR = cgr.getCGR(reac)
        cgrRC = rCGR.get_environment(rCGR.get_center_atoms(), dante=False, deep=2)
        strRC = cgrRC.get_fear(cgrRC.get_morgan())
        if strRC not in dictUnq:
            out.write(cgrRC)
        dictUnq[strRC].append(num)

nameDirect = sys.argv[3]
for i, (k, v) in enumerate(dictUnq.items(), start=1):
    print(k)
    with open('{}/Type_{}.rdf'.format(nameDirect, i), 'w') as fw:
        out = RDFwrite(fw)
        for j in v:
            out.write(reactions[j-1])
"""
