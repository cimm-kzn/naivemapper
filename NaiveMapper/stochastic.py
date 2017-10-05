import math as m
import numpy as np


class Stochastic(object):
    def __init__(self, err=1E-08):
        self.__er = err

    def get_matrix(self, matrix):
        i = 1
        variance = max(max([abs(1-sum(matrix.loc[x])) for x in matrix.index.tolist()]),
                       max([abs(1-sum(matrix[x])) for x in matrix.columns.values.tolist()]))

        def scaling(num):
            if num % 2:
                for j in matrix.index.tolist():
                    matrix.loc[j] /= sum(matrix.loc[j])
                return max([abs(1-sum(matrix[x])) for x in matrix.columns.values.tolist()])
            else:
                for j in matrix.columns.values.tolist():
                    matrix[j] /= sum(matrix[j])
                return max([abs(1-sum(matrix.loc[x])) for x in matrix.index.tolist()])

        while i < 10000 and variance > self.__er:
            variance = scaling(i)
            i += 1

        print(i, '\t', variance)
        for p in matrix.index.tolist():
            for s in matrix.columns.values.tolist():
                # print(matrix.loc[p, s])
                try:
                    matrix.loc[p, s] = - m.log(matrix.loc[p, s])
                except:
                    matrix.loc[p, s] = np.inf

        return matrix
