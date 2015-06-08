__author__ = 'valiaafo'

import numpy as np
from sklearn.naive_bayes import BernoulliNB

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
        BernoulliNB(alpha=1.0, binarize=None, class_prior=None, fit_prior=True)
        log_of_prob = model.predict_log_proba(new_x)
        prob = model.predict_proba(new_x)
        return log_of_prob
        #return prob

    def mapping(self,index_bit,probabilities): #new_y,
        matrix_for_all_reactions = []
        prob_matrix = []
        print(index_bit)
        num = float('inf')
        for i in range(len(index_bit)):
            if index_bit[i][0] != index_bit[i+1][0]:
                matrix_for_all_reactions.append(prob_matrix)
            else:
                for k in range(len(prob_matrix)):
                    for l in range(len(prob_matrix[k])):
                        prob_matrix.append(probabilities[j][0])





