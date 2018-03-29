import pickle
import pandas as pd
from sklearn.naive_bayes import BernoulliNB
from sklearn.neural_network import MLPClassifier

from CGRtools.files.RDFrw import RDFread
from ..bitstringen import Bitstringen
from ..core import getXY, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def retrain_core(**kwargs):
    model = pickle.load(kwargs['model'])  # Загружаем уже существующую модель
    fragger = Fragger(f_type=model['fragment_type'], _min=model['min'], _max=model['max'], _deep=model['deep'],
                      _min2=model['min2'], _max2=model['max2'], _fuzzy=model['fuzzy'])
    bitstring = Bitstringen(model['bitstring'], model['length'], model['fragment_count'])
    pairwise = Pairwise(model['pairs'], model['duplicate'])
    m = model['model']

    # подсчитаем кол-во реакций
    c = 0
    for _ in RDFread(open(kwargs['input'], encoding='cp1251')):
        c += 1

    print("Training set descriptor calculation")
    X_train, Y_train = [], []
    # берем по заданному кол-ву реакции из входящего файла
    for num, reaction in enumerate(worker(RDFread(open(kwargs['input'], encoding='cp1251')), kwargs['debug']), start=1):
        for x, y, _ in getXY(reaction, fragger, pairwise, bitstring, model['chunk']):
            X_train.append(x)
            Y_train.append(y)
            if num % kwargs['batch_chunk'] == 0 or num == c:
                # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)
                m.partial_fit(pd.concat(X_train, ignore_index=True),
                              pd.concat(Y_train, ignore_index=True),
                              classes=pd.Series([False, True]))
                X_train.clear(), Y_train.clear()

    model['model'] = m

    pickle.dump(model, kwargs['model2'])  # записываем нашу обученную модель в файл
