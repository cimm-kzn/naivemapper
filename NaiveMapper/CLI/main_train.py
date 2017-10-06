import pickle

import networkx as nx
import pandas as pd
from sklearn.naive_bayes import BernoulliNB

from CGRtools.files.RDFrw import RDFread
from ..bitstringen import Bitstringen
from ..core import getXY, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def train_core(**kwargs):
    fragger = Fragger(kwargs['min'], kwargs['max'], kwargs['deep'], kwargs['fragment_type'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise = Pairwise(kwargs['pairs'], kwargs['duplicate'])

    # Создаем новую модель Наивного Бейсовского классификатора
    nb = BernoulliNB(alpha=1.0, binarize=None)

    chunk = 200 if kwargs['bitstring'] == 5 else None

    print("Training set descriptor calculation")
    with open(kwargs['input'], encoding='cp1251') as fr:
        for reaction in worker(RDFread(fr)):  # берем по 1 реакции из входящего файла
            for x, y, _ in getXY(reaction, fragger, pairwise, bitstring, chunk):
                # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)
                nb.partial_fit(x, y, classes=pd.Series([False, True]))

    with open(kwargs['model'], 'wb') as train:
        pickle.dump(nb, train)  # записываем нашу обученную модель в файл
