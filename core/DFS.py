from collections import defaultdict
from itertools import chain
from CGRtools.FEAR import FEAR  # from core.Morgan import Morgan

class DFS(object):
    def __init__(self):
        self.__fear = FEAR()

    def getMap(self, s_graph, p_graph, _map):
        new_map = {}
        """
        Создание словаря, где ключом является атом реагента, а значением список атомов в продукте,
        на который может отобразиться данный атом реагента
        """
        p_morgan = self.__fear.get_morgan(p_graph)
        # словарь, ключем которого является нода атом продукта, а занчение его вес(по алгоритму Моргана)
        tmp = defaultdict(list)
        for k, v in p_morgan.items():
            tmp[v].append(k)
            # словарь: ключ - вес, значение - список атомов с одинаковым весом(симметрично-эквивалентные)
        twins = {_map[y]: x for x in tmp.values() for y in x}


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
                #pp_neigh = p_graph.neighbors(new_map[atom])
            except:
                pp_neigh = []

            try:
                return next((x for x in drug if x in pp_neigh and x not in new_map.values()),
                        min(set(drug).difference(new_map.values())))
            except:
                print('Error')

        def walk(trace, inter):
            t = twins[inter]  # Для рассматриваемого атома реагента, находим список сим/экв атомов в продуктах
            p_map = t[0] if len(t) == 1 else morphius([x for x in s_graph.neighbors(inter) if x in new_map], sorted(t))
            # Если список состоит из одного значения, то выбираем его, иначе проводим доп.проверку
            new_map[inter] = p_map  # Новый мап - реагент:продукт
            iterlist = sorted(set(s_graph.neighbors(inter)).difference(trace))
            # Сортируем значения вершин графа субстрата, которые не входили в список ранее рассмотренных атомов
            for i in iterlist:
                if i in trace:
                    continue

                deep = walk(set(chain(trace, [i])), i)
                trace.update(deep)

            return trace

        trace = set(s_graph)
        while trace:
            first = min(trace)
            trace.difference_update(walk({first}, first))

        return {p: s for s, p in new_map.items()}
