import pickle
import pandas as pd
from multiprocess import Pool, Process, Queue
from sklearn.naive_bayes import BernoulliNB
from sklearn.neural_network import MLPClassifier

from CGRtools.files.RDFrw import RDFread
from ..bitstringen import Bitstringen
from ..core import getXY, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def train_core(**kwargs):
    fragger = Fragger(kwargs['min'], kwargs['max'], kwargs['deep'], kwargs['fragment_type'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise = Pairwise(kwargs['pairs'], kwargs['duplicate'])
    model = {k: v for k, v in kwargs.items() if k not in ('input', 'model')}
    # Создаем новую модель Наивного Бейсовского классификатора
    if kwargs['type_model'] is 'nb':
        m = BernoulliNB(alpha=1.0, binarize=None)
    else:
        hls, a, s, alpha = tuple(kwargs['mlp_hls']), kwargs['mlp_a'], kwargs['mlp_s'], kwargs['mlp_alpha']
        print('hidden_layer_sizes:\t{}\nactivation:\t{}\nsolver:\t{}\nalpha:\t{}'.format(hls, a, s, alpha))
        m = MLPClassifier(hidden_layer_sizes=hls, activation=a, solver=s, alpha=alpha)

    # подсчитаем кол-во реакций
    c = 0
    for _ in RDFread(open(kwargs['input'])).read():
        c += 1

    print("Training set descriptor calculation")
    X_train, Y_train = [], []
    # берем по заданному кол-ву реакции из входящего файла
    for num, reaction in enumerate(worker(RDFread(open(kwargs['input'])), kwargs['debug']), start=1):
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

    pickle.dump(model, kwargs['model'])  # записываем нашу обученную модель в файл
