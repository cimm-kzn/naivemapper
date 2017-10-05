import numpy as np

class DFSdb(object):
    def getMap(self, s_graph, p_graph, _map, matrix):
        # res = []
        self.__cost, self.new_map = 0, dict()
        twins = {x: set(m.loc[lambda m_row: abs(m_row - c) < .0000001].index.tolist())
                 for x, m, c in ((x, matrix[x], matrix.loc[y, x]) for y, x in _map.items())}

        def dinamic_bonds(c, s_edge, p_edge, maps):
            '''
            trace = set()
            for s in set(s_edge).intersection(set(maps)):
                if maps[s] not in p_edge:
                    c += 1
                # else:
                    # c += abs(s_edge[s]['s_bond'] - p_edge[maps[s]]['s_bond'])
                trace.add(maps[s])

            se, pe = [s['s_bond'] for s in s_edge.values()], [p['s_bond'] for p in p_edge.values()]
            if se and pe:
                c += abs(np.mean(se) - np.mean(pe))
            else:
                c += np.mean(se)
            # c += abs(sum([s['s_bond'] for s in s_edge.values()]) - sum([p['s_bond'] for p in p_edge.values()]))
            c += len(set(maps.values()).intersection(set(p_edge).symmetric_difference(trace)))
            return c
            '''
            ms, mp = set(maps.keys()) & set(s_edge.keys()), set(maps.values()) & set(p_edge.keys())
            mps = set([s for s, p in maps.items() if p in mp])
            c += len(ms ^ mps)
            c += abs(sum([s_edge[i]['s_bond'] for i in ms]) - sum([p_edge[j]['s_bond'] for j in mp]))
            return c

        def dfs(s_atom, n_map, c):
            for p_atom in sorted(twins[s_atom]-set(n_map.values())):
                c = dinamic_bonds(c, s_graph.adj[s_atom], p_graph.adj[p_atom], n_map)
                if c < self.__cost:
                    n_map[s_atom] = p_atom
                    atom = next((x for x in twins if x > s_atom), False)
                    if atom:
                        dfs(atom, n_map.copy(), c)
                    else:
                        self.new_map = n_map
                        self.__cost = c
                        # res.append(n_map)

        for p, s in _map.items():
            self.new_map[s] = p
            self.__cost = dinamic_bonds(self.__cost, s_graph.adj[s], p_graph.adj[p], self.new_map)
        # print("1st dfs cost={}, map is {}".format(self.__cost, self.new_map))
        dfs(min(twins), dict(), 0)
        # print("Next dfs cost={}, map is {}".format(self.__cost, res))
        return {p: s for s, p in self.new_map.items()}
