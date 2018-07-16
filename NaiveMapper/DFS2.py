import time


class DFSdb(object):
    def getMap(self, s_graph, p_graph, _map, matrix, weights):
        self.__time = time.perf_counter()
        self.__cost, self.new_map = 0, dict()
        # bond_cost = {1: {1: 0, 2: 1, 3: 2, 4: 0}, 2: {1: 1, 2: 0, 3: 1, 4: 0},
        # 3: {1: 2, 2: 1, 3: 0, 4: 1}, 4: {1: 0, 2: 0, 3: 1, 4: 0}}
        twins = {x: set(m.loc[lambda m_row: abs(m_row - c) < .0000001].index.tolist())
                 for x, m, c in ((x, matrix[x], matrix.loc[y, x]) for y, x in _map.items())}

        def dinamic_bonds(c, s_edge, p_edge, s_implicit_h, p_implicit_h, maps):
            ms, mp = set(maps.keys()) & set(s_edge.keys()), set(maps.values()) & set(p_edge.keys())
            mps = set([s for s, p in maps.items() if p in mp])
            c += weights[0]*len(ms ^ mps)
            c += weights[1]*abs(sum([s_edge[i]['s_bond'] for i in ms]) - sum([p_edge[j]['s_bond'] for j in mp]))
            c += weights[2]*abs(s_implicit_h - p_implicit_h)
            return c

        def dfs(s_atom, n_map, c):
            for p_atom in sorted(twins[s_atom]-set(n_map.values())):
                c = dinamic_bonds(c, s_graph.adj[s_atom], p_graph.adj[p_atom],
                                  s_graph.atom_implicit_h(s_atom), p_graph.atom_implicit_h(p_atom), n_map)
                if c < self.__cost:
                    n_map[s_atom] = p_atom
                    atom = next((x for x in twins if x not in n_map), False)
                    t = time.perf_counter() - self.__time
                    if atom and t < 30.0:
                        dfs(atom, n_map.copy(), c)
                    elif t < 30.0 and len(n_map) == len(self.new_map):
                        self.new_map = n_map
                        self.__cost = c

        for p, s in _map.items():
            self.new_map[s] = p
            self.__cost = dinamic_bonds(self.__cost, s_graph.adj[s], p_graph.adj[p],
                                        s_graph.atom_implicit_h(s), p_graph.atom_implicit_h(p), self.new_map)

        dfs(min(twins), dict(), 0)
        return {p: s for s, p in self.new_map.items()}
