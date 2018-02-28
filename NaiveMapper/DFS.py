from collections import defaultdict
from itertools import chain
# from CGRtools.strings import get_morgan


def equivalent(morgan_dict):
    """ словарь, ключем которого является нода атома, а занчение его вес(по алгоритму Моргана)
    превратить в словарь: ключ - нода атома, значение - список атомов с одинаковым весом(симметрично-эквивалентные)
    """
    tmp = defaultdict(list)
    for k, v in morgan_dict.items():
        tmp[v].append(k)
        # словарь: ключ - вес, значение - список атомов с одинаковым весом(симметрично-эквивалентные)
    return {y: x for x in tmp.values() for y in x}


def get_map_dfs(s_graph, p_graph, _map, *_):
    """
    s_grup = equivalent(self.__fear.get_morgan(s_graph))
    # {атом_реагента: [атомы_реагента, сим.\экв. с рассмтриваемым по Моргану)]}
    twins = {x: m.loc[lambda m_row: abs(m_row - c) < .0000001].index.tolist()
             for x, m, c in ((x, matrix[x], matrix.loc[y, x]) for y, x in _map.items())}
    # {атом_реагента:[атомы_продукта, на который может отобразиться данный атом реагента (по значениям Likelihood)]}
    """

    def morphius(atoms, drug):
        """
        Если есть уже смапированные атомы соседи, то у них находятся общий сосед в атомах продуктах.
        Если нет, то создается пустой список. Выбор среди группы сим/экв атомов: минимальное значение ноды,
        которая входит в созданный список и не входило в ранее смапированный список атомов
        """
        try:
            pp = [set(p_graph.neighbors(new_map[x])) for x in atoms]
            if len(pp) > 1:
                for i in range(1, len(pp)):
                    pp[0].intersection_update(pp[i])
            pp_neigh = pp[0]
        except:
            pp_neigh = []

        try:
            return next((x for x in drug if x in pp_neigh and x not in new_map.values()),
                        min(set(drug).difference(new_map.values())))
        except:
            print('Error')

    def walk(trace, inter):
        t1 = p_grup[inter]
        # t2 = s_grup[inter]
        # Для рассматриваемого атома реагента, находим список сим/экв атомов в продуктах и реагентах
        p_map = morphius([x for x in s_graph.neighbors(inter) if x in new_map], sorted(t1))
        #                 sorted(twins[inter]) if len(t1) == 1 and len(t2) != 1 else sorted(t1))
        # Если нет сим/экв атомов в продуктах, но они имеются в реагенте, то рассматриваем атомы с равным Likelihood
        new_map[inter] = p_map  # Новый мап - реагент:продукт
        iterlist = sorted(set(s_graph.neighbors(inter)).difference(trace))
        # Сортируем значения вершин графа субстрата, которые не входили в список ранее рассмотренных атомов
        for i in iterlist:
            if i in trace:
                continue

            deep = walk(set(chain(trace, [i])), i)
            trace.update(deep)

        return trace

    new_map = {}
    p_grup = {_map[y]: x for y, x in equivalent(p_graph.get_morgan()).items()}
    # {атом_реагента: [атомы_продукта, на который может отобразиться данный атом реагента (по Моргану)]}

    trace = set(s_graph)
    while trace:
        first = min(trace)
        trace.difference_update(walk({first}, first))

    return {p: s for s, p in new_map.items()}
