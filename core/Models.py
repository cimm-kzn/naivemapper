__author__ = 'valiaafo'

import numpy as np
from sklearn.naive_bayes import BernoulliNB
from sklearn.utils.tests.test_linear_assignment import _hungarian

class Models(object):
    def __init__(self):
        pass


    def learning(self,x,y,**parameters):
        X = x
        Y = y
        model = BernoulliNB()
        model.fit(X, Y)
        print(model)
        return model


    def predict(self,model,new_x):
        # BernoulliNB(alpha=1.0, binarize=None, class_prior=None, fit_prior=True)
        log_of_prob = model.predict_log_proba(new_x)
        print(log_of_prob)
        prob = model.predict_proba(new_x)
        return log_of_prob
        #return prob

    def mapping(self,index_bit,probabilities,quantity_a):
        print(index_bit)
        num = float('inf')
        prob_matrix = np.zeros((quantity_a,quantity_a))
        print(prob_matrix)
        for i in range(len(probabilities)):
            pr_i = probabilities[i][0]*(-1)
            k = index_bit[i][1]
            l = index_bit[i][2]
            prob_matrix[k][l] = pr_i
        for m in range(len(prob_matrix)):
            for n in range(len(prob_matrix[m])):
                if prob_matrix[m][n] == 0:
                    prob_matrix[m][n] = num
        print(prob_matrix)
        indexes = _hungarian(prob_matrix)
        total_cost = 0
        for s, p in indexes:
            x = prob_matrix[s, p]
            total_cost += x
        print(indexes)
        print(total_cost)







