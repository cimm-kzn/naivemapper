from typing import Tuple, Iterable
from itertools import product
from CGRtools.strings import get_morgan
import networkx as nx
import pandas as pd
import copy


class Pairwise(object):
    def __init__(self, p_type=0, d_type=True):
        self.__pairs_type = self.equivalent if p_type else self.simple
        if p_type:
            self.__duplicate_type = self.__doesnt if d_type else self.__has

    def get(self, sub_graph: nx.Graph, prod_graph: nx.Graph) -> (Iterable[Tuple[int, int]], pd.Series):
        """
        Возвращает номера пар одинаковых атомов реагента и продукта (возможны разные принципы создания пар).
        Т.е. для разных типов атомов, N и C, пары не целесообразно создавать
        """
        return self.__pairs_type(sub_graph, prod_graph)

    @staticmethod
    def simple(sub_graph, prod_graph):
        pairs = []
        state = []
        for (s_atom, s_prop), (p_atom, p_prop) in product(sub_graph.nodes(data=True), prod_graph.nodes(data=True)):
            if s_prop['element'] == p_prop['element']:
                pairs.append((s_atom, p_atom))
                state.append(s_atom == p_atom)
        return pairs, pd.Series(state)

    def equivalent(self, sub_graph, prod_graph):
        sub_m, prod_m = get_morgan(sub_graph), get_morgan(prod_graph)
        s_grup, p_grup = {x: [] for x in set(sub_m.values())}, {y: [] for y in set(prod_m.values())}
        for (k1, v1), (k2, v2) in zip(sub_m.items(), prod_m.items()):
            s_grup[v1].append(k1)
            p_grup[v2].append(k2)

        grup = list(product(s_grup.keys(), p_grup.keys()))
        g = copy.copy(grup)
        c = []
        for i in g:
            list1 = s_grup[i[0]]
            list2 = p_grup[i[1]]

            if sub_graph.node[list1[0]]['element'] == prod_graph.node[list2[0]]['element']:
                c.append(len(list1)-len(set(list1)-set(list2)) != 0)
            else:
                grup.remove(i)

        return self.__duplicate_type(grup, c, s_grup, p_grup)

    @staticmethod
    def __doesnt(grup, c, s_grup, p_grup):
        pairs = []
        for i in grup:
            pairs.append((s_grup[i[0]][0], p_grup[i[1]][0]))

        return pairs, pd.Series(c)

    @staticmethod
    def __doesFalse(grup, c, s_grup, p_grup):  # Not used
        y = []
        pairs = []
        for j, i in enumerate(grup):
            if c[j]:
                p = list(product(s_grup[i[0]], p_grup[i[1]]))
                for k in range(len(p)):
                    y.append(c[j])
                    pairs.append(p[k])
            else:
                y.append(c[j])
                pairs.append((s_grup[i[0]][0], p_grup[i[1]][0]))

        return pairs, pd.Series(y)

    @staticmethod
    def __has(grup, c, s_grup, p_grup):
        y = []
        pairs = []
        for j, i in enumerate(grup):
            p = list(product(s_grup[i[0]], p_grup[i[1]]))
            for k in range(len(p)):
                y.append(c[j])
                pairs.append(p[k])

        return pairs, pd.Series(y)
