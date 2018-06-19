from keras.utils import to_categorical
import pickle
import pandas as pd
from sklearn.naive_bayes import BernoulliNB
from sklearn.neural_network import MLPClassifier
from ..keras_class.keras_mlp import Keras_MLP
# from ..keras_mlp_class.keras_mlp import Keras_MLP
from keras import backend as kb

from CGRtools.files.RDFrw import RDFread
from ..bitstringen import Bitstringen
from ..core import getXY, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


kb.set_session(kb.tf.Session(config=kb.tf.ConfigProto(inter_op_parallelism_threads=2, intra_op_parallelism_threads=2)))

def train_core(**kwargs):
    fragger = Fragger(f_type=kwargs['fragment_type'], _min=kwargs['min'], _max=kwargs['max'], _deep=kwargs['deep'],
                      _min2=kwargs['min2'], _max2=kwargs['max2'], _fuzzy=kwargs['fuzzy'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise = Pairwise(kwargs['pairs'], kwargs['duplicate'])
    model = {k: v for k, v in kwargs.items() if k not in ('input', 'model')}
    # Создаем новую модель Наивного Бейсовского классификатора
    if kwargs['type_model'] == 'nb':
        m = BernoulliNB(alpha=1.0, binarize=None)
    elif kwargs['type_model'] == 'mlp':
        hls, a, s, alpha = tuple(kwargs['mlp_hls']), kwargs['mlp_a'], kwargs['mlp_s'], kwargs['mlp_alpha']
        bs, lr, mi, es = kwargs['mlp_bs'], kwargs['mlp_lr'], kwargs['mlp_mi'], kwargs['mlp_es']
        m = MLPClassifier(hidden_layer_sizes=hls, activation=a, solver=s, alpha=alpha,
                          batch_size=bs, learning_rate=lr, max_iter=mi, early_stopping=es)
    else:
        bs, hls, alpha, es = kwargs['mlp_bs'], tuple(kwargs['mlp_hls']), kwargs['mlp_alpha'], kwargs['mlp_es']
        if kwargs['keras_dropout']: do = kwargs['keras_dropout']
        else: do = "Auto"

        classifier = Keras_MLP(task="classification", layer_sizes=hls, activations='relu', dropout=do, alpha=alpha,
                               batch_size=bs, learning_rate_init=0.001, epochs=1, shuffle=True, loss_function="mse",
                               metrics=['mse'], verbose=kwargs['debug'], early_stopping=es, optimizer_name="adam",
                               lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
        tb = [1, 3, 3, 2, 2, 1]
        m = classifier.create_model(kwargs['length']*tb[kwargs['bitstring']], 2)
        # print(m.layer_sizes)

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
            # if num == 1 and kwargs['type_model'] == 'keras':
            #     m = classifier.create_model(x.values.shape[1], 1)
            if num % kwargs['batch_chunk'] == 0 or num == c:
                # обучаем нашу модель на основании битовыx строк дескрипторов(Х) и строк значений ААО(Y)

                if kwargs['type_model'] == 'keras':
                    """
                    x_np = pd.concat(X_train, ignore_index=True).values
                    y_np = to_categorical(1*pd.DataFrame(pd.concat(Y_train, ignore_index=True)).values)
                    m.partial_fit(x_np, y_np, kwargs['model_filename'])
                    """
                    m.fit(pd.concat(X_train, ignore_index=True).values,
                          to_categorical(1 * pd.DataFrame(pd.concat(Y_train, ignore_index=True)).values),
                          batch_size=classifier.batch_size, epochs=classifier.epochs, verbose=classifier.verbose,
                          callbacks=classifier.used_callbacks, shuffle=classifier.shuffle)
                else:
                    m.partial_fit(pd.concat(X_train, ignore_index=True),
                                  pd.concat(Y_train, ignore_index=True),
                                  classes=pd.Series([False, True]))

                X_train.clear(), Y_train.clear()
    if kwargs['type_model'] == 'keras':
        model['model_filename'] = kwargs['model_filename']
        m.save(kwargs['model_filename'])
        pickle.dump(model, kwargs['model'])  # записываем настроенные параметры модели
    else:
        model['model'] = m
        pickle.dump(model, kwargs['model'])  # записываем нашу обученную модель в файл
