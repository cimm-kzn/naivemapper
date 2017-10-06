#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
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
from collections import Counter
from CGRtools.strings import get_cgr_string
from itertools import tee
from typing import Dict, Tuple
import networkx as nx


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ...
    copypaste from off-doc
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Fragger(object):
    def __init__(self, _min: int = 1, _max: int = 10, _deep: int = 3, f_type: int=2):
        self.__min = _min
        self.__max = _max
        self.__deep = _deep
        self.__type = f_type
        self.__fragment = self.__sequences if self.__type == 0 else self.__augmented if self.__type == 1 \
            else self.__aug_and_seq

    def get(self, data: nx.Graph) -> Dict[int, Dict[Tuple[Tuple], int]]:
        return self.__fragment(data)

    def __sequences(self, data):
        """
        создает словарь атомов и словарей фрагментов (цепочки атомов и связей) и их количеств.
        """
        def path_namer(path):
            frag = []
            for x in path:
                node_info = data.node[x]
                # frag.append('#%s%s' % (node_info['element'],node_info['s_charge']))
                # убрали заряды из-за того, что общие фрагменты единочной длины исчезали (например '#N0'!='#N3')
                frag.append('#%s' % (node_info['element']))
            for x, y in pairwise(path):
                edge_info = data[x][y]
                frag.append(str(edge_info['s_bond']))

            return '.'.join(frag)

        paths = nx.all_pairs_shortest_path(data, cutoff=self.__max - 1)
        return {x: dict(Counter(path_namer(z) for z in y.values() if len(z) >= self.__min)) for x, y in paths}

    def __augmented(self, data):
        """
        создает словарь атомов и словарей фрагментов (атомы с окружением/augmented atoms).
        """
        return {x: dict.fromkeys({'%d^%s' % (n, get_cgr_string(y, y.get_morgan()))
                                  for n, y in enumerate(data.get_environment([x], dante=True, deep=self.__deep))
                                  if (n and self.__type == 2) or self.__type in [0, 1]}, 1)
                for x in data.nodes()}

    def __aug_and_seq(self, data):
        """Создает список фрагментов состоящих из augmented(1:n) и sequences(0:3)"""
        aug = self.__augmented(data)
        return {x: dict(y, **aug[x]) for x, y in self.__sequences(data).items()}
