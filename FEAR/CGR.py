#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Find Errors in Automapped Reactions).
#
# FEAR is free software; you can redistribute it and/or modify
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
import copy
import numpy
from mendel import replace
from maprules import rules


def main():
    print("This file is part of FEAR.")
    return 0


class CGR(object):
    def __init__(self):
        self.__replace = replace()
        self.__rulesC = rules('rules_c.dat')
        self.__rulesE = rules('rules_e.dat')

    def firstcgr(self, data, rid):
        '''
        генерит матрицы продуктов и реагентов.
        сортирует их и дополняет атомами приведенными только в продуктах или реагентах.
        :param data:  словарь описывающий всю реакцию
        :return:
        '''
        self.__rid = rid
        self.__data = {}
        #self.__mols = dict(substrats=[], products=[])
        self.__data['substrats'] = self.getmatrix(data['substrats'])
        self.__data['products'] = self.getmatrix(data['products'])
        self.__sortmatrix()
        self.__mapsrepl = {}
        #try:
        self.__scan()
        #data['diff'] = self.__diff()
        rc = self.__reactcenter()
        #rv = numpy.absolute(self.__data['products']['bondmatrix'] - self.__data['substrats']['bondmatrix']).sum()
        for i, j in enumerate(rc[1]):
            data['meta']['!reaction_center_hash_%d' % i] = j
        return rc[0]

    def __diff(self):
        diff = []
        for x, y in self.__data['diff'].items():
            for w, z in y.items():
                if (w + 1, x + 1, z) not in diff:
                    diff += [(x + 1, w + 1, z)]
        return diff

    def getmatrix(self, data):
        '''
        объединение матриц связности группы молекул
        '''
        atomlist = []
        #объединим атомы в 1 список. и создадим пустую матрицу.
        for i in data:
            #self.__mols[role].append(range(len(atomlist), len(atomlist) + len(i['atomlist'])))
            maps = [x['map'] for x in atomlist]
            collide = set([x['map'] for x in i['atomlist']]).intersection(maps)
            if collide and collide != {0}:
                print set([x['map'] for x in i['atomlist']]).intersection(maps), self.__rid
                rep = max(maps)
                for x in i['atomlist']:
                    x['map'] += rep
            atomlist += i['atomlist']
        length = len(atomlist)
        self.__length = length
        tempMatrix = {'bondmatrix': numpy.zeros((length, length), dtype=float), 'atomlist': atomlist}
        step = 0
        for i in data:
            matrixlength = len(i['atomlist'])
            tempMatrix['bondmatrix'][step:step + matrixlength, step:step + matrixlength] = i['bondmatrix']
            step += matrixlength
        for i in range(length):#zip(*tempMatrix.values()):
            j = list(tempMatrix['bondmatrix'][i])
            j.pop(i)
            if 1.5 in j:
                t = 4
            elif 3 in j:
                t = 3
            elif 2 in j:
                t = 2
            else:
                t = 1
            atomlist[i]['type'] = t

        return tempMatrix

    def __sortmatrix(self):
        '''
        атомы отсортировываются по порядку. для атомов представленных только в продуктах или реагентах достраивается
        симметричная матрица с разрывом в месте стыка с основной молекулой.
        '''
        maps = {'substrats': [x['map'] for x in self.__data['substrats']['atomlist']],
                'products': [x['map'] for x in self.__data['products']['atomlist']]}

        length = max((max(maps['products']), max(maps['substrats'])))

        for j in ('substrats', 'products'):
            if 0 in maps[j]:
                nmaps = []
                for x in maps[j]:
                    if x:
                        nmaps.append(x)
                    else:
                        length += 1
                        nmaps.append(length)
                maps[j] = nmaps

        lose = sorted(list(set(range(1, length + 1)).difference(maps['products']).difference(maps['substrats'])),
                      reverse=True)
        if lose:
            for k in maps.keys():
                for i in lose:
                    newmaps = []
                    for j in maps[k]:
                        newmaps.append(j if j < i else j - 1)
                    maps[k] = newmaps
            length -= len(lose)

        self.__length = length

        for j in ('substrats', 'products'):
            tmp = numpy.zeros((length, length), dtype=float)

            shape = self.__data[j]['bondmatrix'].shape[0]
            for i in xrange(len(maps[j])):
                tmp[maps[j][i] - 1, 0:shape] = self.__data[j]['bondmatrix'][i, :]
            self.__data[j]['bondmatrix'] = numpy.zeros((length, length), dtype=float)
            for i in xrange(len(maps[j])):
                self.__data[j]['bondmatrix'][:, maps[j][i] - 1] = tmp[:, i]

            atomlist = [0] * self.__length
            for i, k in zip(self.__data[j]['atomlist'], maps[j]):
                i['map'] = k
                atomlist[k - 1] = i
            self.__data[j]['atomlist'] = atomlist
        #if 0 in self.__data['substrats']['atomlist'] or 0 in self.__data['products']['atomlist']:
        #    print "ACHTUNG"

        if 0 in self.__data['substrats']['atomlist']:
            '''
            метод восстановления используется только для недостающих атомов реагентов.
            информация о факте отвалившегося фрагмента в продуктах самодостаточна для корректного ААО или
            восстановления структуры уходящей группы.
            '''
            self.__moltorepare = True
            self.__searchlost('products', 'substrats')
        if 0 in self.__data['products']['atomlist']:
            self.__searchlost('substrats', 'products')

    def __searchlost(self, s, t):
        '''
        внедрение в матрицу продуктов/субстратов атомов из субстратов/продуктов
        '''
        lostlist = []
        for i in xrange(self.__length):
            if self.__data[t]['atomlist'][i] == 0:
                self.__data[t]['atomlist'][i] = self.__data[s]['atomlist'][i]
                lostlist.append(i)
        for i in lostlist:
            for j in lostlist:
                self.__data[t]['bondmatrix'][i, j] = self.__data[s]['bondmatrix'][i, j]

    def __reactcenter(self):
        atoms = []
        report = []
        newstruct = 0

        def getsmarts(sub, level):
            resout = []
            while sub:
                self.__nextnumb = self.__numb()
                res = self.__getSmiles(sub[:1], sub[0], level=level)
                resout.append(''.join(res[1]))
                sub = list(set(sub) - res[0])
            return resout

        rv = numpy.absolute(self.__data['products']['bondmatrix'] - self.__data['substrats']['bondmatrix'])
        for i in xrange(self.__length):
            if i not in atoms:
                sub = list(self.__reactcenteratoms([i], i))
                atoms.extend(sub)

                if len(sub) > 1:
                    srv = rv[numpy.ix_(sub, sub)].sum()
                    hashes = []
                    for j in range(2, -1, -1):
                        mhash = self.__getMorgan(sub, level=j)
                        mhashkey = tuple(sorted(mhash.values()))
                        ruledataC = self.__rulesC.get(mhashkey, [0])
                        ruledataE = self.__rulesE.get(mhashkey, [0])
                        if ruledataC[0]:
                            hashes = [ruledataC[1]]
                            newhash = 0
                            break
                        elif ruledataE[0]:
                            hashes = [ruledataC[1]]
                            newhash = 1
                            break
                        else:
                            newhash = 1
                            self.__recsub = sub
                            self.__smartstarget = 'substrats'
                            self.__nextmap = self.__numb()
                            smirks = '.'.join(getsmarts(sub, j))
                            self.__smartstarget = 'products'
                            smirks += '>>' + '.'.join(getsmarts(sub, j))
                            hashes.append("new: %s'%.1f + %s" % (','.join(str(x) for x in mhashkey), srv, smirks))
                    report.extend(hashes)
                    newstruct += newhash
        return (True, report) if not newstruct else (False, report)

    __tosmiles = {1: '-', 2: '=', 3: '#', 1.5: ':'}

    def __numb(self):
        i = 1
        while True:
            yield i
            i += 1

    def __getSmiles(self, trace, inter, level=None):
        strace = set(trace)
        rule = {'substrats': lambda x: x < 80, 'products': lambda x: x % 10 != 8}
        iterlist = set([x for x, y in self.__data['diff'][inter].items() if rule[self.__smartstarget](y)]).intersection(self.__recsub).difference(trace[-2:])
        #print iterlist, inter, self.__data['diff'][inter]
        if level == 1:
            symbol = '#%d' % self.__replace[self.__data[self.__smartstarget]['atomlist'][inter]['element']]
        elif level == 2: # с типами атомов
            symbol = '#%d%s' % (self.__replace[self.__data[self.__smartstarget]['atomlist'][inter]['element']], self.__types[self.__data[self.__smartstarget]['atomlist'][inter]['type']])
        else:
            symbol = '*'
        if self.__smartstarget == 'substrats':
            self.__mapsrepl[inter + 1] = self.__nextmap.next()
        smi = ['[%s:%d]' % (symbol, self.__mapsrepl[inter + 1])]
        #print smi
        concat = []
        stoplist = []
        for i in iterlist:
            if i in strace:
                if i not in stoplist:
                    # костыль для циклов. чтоб не было 2х проходов.
                    cyc = self.__nextnumb.next()
                    concat += [(i, cyc, inter)]
                    #print concat
                    smi[0] += '%s%d' % (self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]], cyc)
                    #print smi, '!'
                continue
            deep = self.__getSmiles(copy.copy(trace + [i]), i, level)
            strace.update(deep[0])
            #print strace
            #print deep
            if deep[2]:
                #print inter, i, '1'
                concat += deep[2]
                #print deep[2]
                for j in deep[2]:
                    if j[0] == inter:
                        #print '2', inter
                        stoplist += [j[2]]
                        #print stoplist
                        smi[0] += '%s%d' % (self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]], j[1])
            smi += ['(%s' % self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]]] + deep[1] + [')']

            #strace.update(self.__reactcenteratoms(copy.copy(trace + [i]), i))
        return strace, smi, concat

    def __reactcenteratoms(self, trace, inter):
        '''
        поиск всех атомов подструктур-реакционных центров
        '''
        strace = {inter}
        iterlist = set(self.__data['diff'][inter]).difference(trace)

        for i in iterlist:
            if i in strace:
                # костыль для циклов. чтоб не было 2х проходов.
                continue
            if self.__data['diff'][inter][i] > 10:
                strace.update(self.__reactcenteratoms(copy.copy(trace + [i]), i))
        return strace

    __types = {4: 'A', 3: 'T', 2: 'D', 1: 'S'}

    def __getMorgan(self, sub, level=None):
        '''
        модифицированный алгоритм моргана используется для хэширования подструктур.
        изначально все атомы имеют веса согласно сорту атомов.
        итеративно значения весов атомов умножатся на кратности образуемых связей с последующим суммированием полученных
        значений к текущему весу.

        после возвращается хэш подструктуры.
        '''
        if level == 1:
            weights = {x: self.__replace[self.__data['products']['atomlist'][x]['element']] * 1000 for x in sub}
        elif level == 2:
            weights = {x: self.__replace[self.__data['products']['atomlist'][x]['element']] * 1000 + self.__data['products']['atomlist'][x]['type'] * 1000000 for x in sub}
        else:
            weights = dict.fromkeys(sub, 0)
        for i in sub:
            l = sum([x for x in self.__data['diff'][i].values() if x > 10]) + self.__data['diff'][i][i]
            weights[i] += l
        return weights

    __cgr = {0: {0: 0, 1: 81, 2: 82, 3: 83, 4: 84, 5: 85, 6: 86, 7: 87, 8: 88},
             1: {0: 18, 1: 1, 2: 12, 3: 13, 4: 14, 5: 15, 6: 16, 7: 17, 8: 18},
             2: {0: 28, 1: 21, 2: 2, 3: 23, 4: 24, 5: 25, 6: 26, 7: 27, 8: 28},
             3: {0: 38, 1: 31, 2: 32, 3: 3, 4: 34, 5: 35, 6: 36, 7: 37, 8: 38},
             4: {0: 48, 1: 41, 2: 42, 3: 43, 4: 4, 5: 45, 6: 46, 7: 47, 8: 48},
             5: {0: 58, 1: 51, 2: 52, 3: 53, 4: 54, 5: 5, 6: 56, 7: 57, 8: 58},
             6: {0: 68, 1: 61, 2: 62, 3: 63, 4: 64, 5: 65, 6: 6, 7: 67, 8: 68},
             7: {0: 78, 1: 71, 2: 72, 3: 73, 4: 74, 5: 75, 6: 76, 7: 7, 8: 78},
             8: {0: 80, 1: 81, 2: 82, 3: 83, 4: 84, 5: 85, 6: 86, 7: 87, 8: 8}}

    def __scan(self):
        self.__data['diff'] = {x: {} for x in range(self.__length)}

        for i in xrange(self.__length):
            s = self.__data['substrats']['bondmatrix'][i, i]
            p = self.__data['products']['bondmatrix'][i, i]
            self.__data['diff'][i][i] = int(p - s)
            for j in xrange(i + 1, self.__length):
                s = self.__data['substrats']['bondmatrix'][i, j]
                p = self.__data['products']['bondmatrix'][i, j]
                if s != 0 or p != 0:
                    diff = self.__cgr[s if s != 1.5 else 4][p if p != 1.5 else 4]
                    self.__data['diff'][i][j] = self.__data['diff'][j][i] = diff

    def __accept(self, numb, length):
        '''
        функция остановки для алгоритма моргана. останавливает поиск после первого простоя улучшения.
        '''
        numl = len(set(numb.values()))
        if numl == length:
            return False
        if self.__buffer == numl:
            return False
        self.__buffer = numl
        return True


if __name__ == '__main__':
    main()
