#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of naivemapper.
#
# naivemapper is free software; you can redistribute it and/or modify
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


def main():
    print("This file is part of naivemapper.")
    return 0


class Mapper(object):
    def __init__(self):
        self.__data = {}

    def get(self, data):
        '''
        генерит матрицы продуктов и реагентов.
        сортирует их и дополняет атомами приведенными только в продуктах или реагентах.
        :param data:  словарь описывающий всю реакцию
        :return:
        '''
        self.__data['substrats'] = self.__getmatrix(data['substrats'])
        self.__data['products'] = self.__getmatrix(data['products'])
        print(self.__data)
        return True

    def __getmatrix(self, data):
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


if __name__ == '__main__':
    main()
