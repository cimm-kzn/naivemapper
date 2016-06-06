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
from itertools import tee
import networkx as nx
from typing import Dict, Tuple
from collections import Counter


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ...
    copypaste from off-doc
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Fragger(object):
    def __init__(self, _min: int = 1, _max: int = 10):
        self.__min = _min
        self.__max = _max

    def get(self, data: nx.Graph) -> Dict[int, Dict[Tuple[Tuple], int]]:
        """
        создает словарь атомов и словарей фрагментов и их количеств.
        """
        def path_namer(path):
            frag = []
            for x in path:
                node_info = data.node[x]
                frag.append((node_info['element'],node_info['s_charge']))
            for x, y in pairwise(path):
                edge_info = data[x][y]
                frag.append((edge_info['s_bond'],))

            return tuple(frag)

        paths = nx.all_pairs_shortest_path(data, cutoff=self.__max - 1)
        return {x: dict(Counter(path_namer(z) for z in y.values() if len(z) >= self.__min)) for x, y in paths.items()}
