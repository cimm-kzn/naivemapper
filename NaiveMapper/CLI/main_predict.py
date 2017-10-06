import pickle

import networkx as nx

from CGRtools.files.RDFrw import RDFread, RDFwrite
from ..bitstringen import Bitstringen
from ..core import getXY, mapping, truth, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def predict_core(**kwargs):
    fragger = Fragger(kwargs['min'], kwargs['max'], kwargs['deep'], kwargs['fragment_type'])
    bitstring = Bitstringen(kwargs['bitstring'], kwargs['length'], kwargs['fragment_count'])
    pairwise = Pairwise(0, 0)  # при предсказании, не применяем алгоритма Моргана для генерации пар сим./экв. атомов

    # Загружаем уже существующую модель в Наивный Бейсовский классификатор
    with open(kwargs['model'], 'rb') as train:
        nb = pickle.load(train)

    chunk = 200 if kwargs['bitstring'] == 5 else None

    print("Testing set descriptor calculation")
    with open(kwargs['input'], encoding='cp1251') as fr, open(kwargs['output']) as fw:
        outputdata = RDFwrite(fw)
        for reaction in worker(RDFread(fr)):  # берем по 1 реакции из входящего файла
            y, pairs = [], []
            for x, _, drop_pairs in getXY(reaction):
                '''на основании сгенерированного набора битовых строк дескрипторов
                из обученой модели выводятся значения лагорифмов вероятностей проецирования (отображения) атомов'''
                pairs.extend(drop_pairs)
                y.extend([yy for yy in nb.predict_log_proba(x)])

            _map = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
            tmp = []
            for graph in reaction['products']:
                tmp.append(nx.relabel_nodes(graph, _map, copy=True))
                # на основании обученной модели перемаппливаются атомы продукта
            reaction['products'] = tmp
            # reaction['meta']['Likelihood'] = opt_lh  # доп.графа в исх.файле, со ср.знач.вероятностей отображения

            '''for s, p_lh in lh.items():
                lh_name = str(s) + '_ATOM_PROBABILITY_MAP'
                reaction['meta'][lh_name] = p_lh
                # доп.графы в исх.файле, с вероятностями отображения данного атома реагента на все атомы продукта
            '''
            outputdata.write(reaction)  # запись реакции в исходящий файл
    _, _ = truth(kwargs['input'], kwargs['output'], 0, 0, [])  # Проверка соответствия

