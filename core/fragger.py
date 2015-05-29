#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
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


def main():
    print("This file is part of naivemapper.")
    return 0


class Fragger(object):
    def __init__(self, **kwargs):
        self.__min = kwargs['min']
        self.__max = kwargs['max']

    def get(self, data):
        '''
        генерит матрицы продуктов и реагентов.
        сортирует их и дополняет атомами приведенными только в продуктах или реагентах.
        :param data:  словарь описывающий всю реакцию
        :return:
        '''
        result = {'substrats': self.__dosearch(self.__getmatrix(data['substrats'])),
                  'products': self.__dosearch(self.__getmatrix(data['products']))}

        return result

    @staticmethod
    def __getmatrix(data):
        '''
        объединение матриц связности группы молекул
        '''
        repl = []
        connections = {}

        for i in data:
            lrepl = len(repl)
            for x, y in enumerate(i['atomlist']):
                key = (x + lrepl, y['element'], y['izotop'])
                connections[key] = {}
                repl.append(key)


            for x in i['bondlist']:
                connections[repl[x[0] + lrepl]][repl[x[1] + lrepl]] = connections[repl[x[1] + lrepl]][
                    repl[x[0] + lrepl]] = x[2]
        return connections

    def __dosearch(self, connections):
        self.__connections = connections
        total = {}
        for i in self.__connections:
            self.__fragments = defaultdict(int)
            self.__lcs([i], i)
            total[i] = self.__fragments
        return total

    def __lcs(self, trace, inter):
        iterlist = list(set(self.__connections[inter]) - set(trace))
        if len(trace) // 2 + 1 < self.__max:
            for i in iterlist:
                self.__lcs(trace + [self.__connections[trace[-1]][i], i], i)
        if len(trace) // 2 + 1 >= self.__min:  # or not iterlist and len(self.__connections[trace[0]]) < 2 and self.__small:
            atoms = [tuple(x[1:]) for x in trace[::2]]
            frag = [item for sublist in zip(atoms, trace[1::2]) for item in sublist] + atoms[-1:]
            self.__fragments[tuple(frag)] += 1

if __name__ == '__main__':
    main()
