import pandas as pd
import pickle
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
    model = {k: v for k, v in kwargs.items() if k not in ('input', 'model')}
    # Создаем новую модель Наивного Бейсовского классификатора
    nb = BernoulliNB(alpha=1.0, binarize=None)

    print("Training set descriptor calculation")
    for reaction in worker(RDFread(kwargs['input']), kwargs['debug']):  # берем по 1 реакции из входящего файла
        for x, y, _ in getXY(reaction, fragger, pairwise, bitstring, kwargs['chunk']):
            # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)
            nb.partial_fit(x, y, classes=pd.Series([False, True]))
    model['model'] = nb

    pickle.dump(model, kwargs['model'])  # записываем нашу обученную модель в файл
