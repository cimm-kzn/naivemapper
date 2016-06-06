from typing import Tuple, Iterable
import networkx as nx
import pandas as pd
from itertools import product as pr_iter


class Wedding(object):
    def __init__(self, p_type=0):
        self.__pairs_type = self.dumb if p_type == 0 else None

    def get(self, sub_graph: nx.Graph, prod_graph: nx.Graph) -> Tuple[Iterable[Tuple[int, int]], pd.Series]:
        """
        Возвращает номера пар одинаковых атомов реагента и продукта (возможны разные принципы создания пар).
        Т.е. для разных типов атомов, N и C, пары не целесообразно создавать
        """
        return self.__pairs_type(sub_graph, prod_graph)

    def dumb(self, sub_graph, prod_graph):
        pairs = []
        state = []
        for (s_atom, s_prop), (p_atom, p_prop) in pr_iter(sub_graph.nodes(data=True), prod_graph.nodes(data=True)):
            if s_prop['element'] == p_prop['element']:
                pairs.append((s_atom, p_atom))
                state.append(s_atom == p_atom)
        return pairs, pd.Series(state)
