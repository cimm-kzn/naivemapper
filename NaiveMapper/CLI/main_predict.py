import networkx as nx
import pickle

from CGRtools.files.RDFrw import RDFread, RDFwrite
from ..bitstringen import Bitstringen
from ..core import getXY, mapping, truth, worker
from ..fragger import Fragger
from ..pairwise import Pairwise


def remap(graphs, maps):
    tmp = []
    for graph in graphs:
        tmp.append(graph.remap(maps, copy=True))
    return tmp


def predict_core(**kwargs):
    # Загружаем уже существующую модель в Наивный Бейсовский классификатор
    model = pickle.load(kwargs['model'])
    fragger = Fragger(model['min'], model['max'], model['deep'], model['fragment_type'])
    bitstring = Bitstringen(model['bitstring'], model['length'], model['fragment_count'])
    # При предсказании, НИКОГДА не применяем алгоритма Моргана для генерации пар сим./экв. атомов
    pairwise = Pairwise(0, False)

    print("Testing set descriptor calculation")
    with open(kwargs['input'], encoding='cp1251') as fr, open(kwargs['output'], 'w') as fw:
        outputdata = RDFwrite(fw)
        for reaction in worker(RDFread(fr), kwargs['debug']):
            p_graph = nx.union_all(reaction['products'])
            reaction['products'] = remap(reaction['products'], {k: k + max(p_graph.nodes()) for k in p_graph.nodes()})
            y, pairs = [], []
            for x, _, drop_pairs in getXY(reaction, fragger, pairwise, bitstring, model['chunk']):
                '''
                На основании сгенерированного набора битовых строк дескрипторов из обученой модели 
                выводятся значения логорифмов вероятностей проецирования (отображения) атомов
                '''
                pairs.extend(drop_pairs)
                # y.extend([yy for yy in nb.predict_log_proba(x)])
                y.extend(model['model'].predict_log_proba(x))

            _map = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
            # на основании обученной модели перемаппливаются атомы продукта
            reaction['products'] = remap(reaction['products'], _map)

            '''
            # reaction['meta']['Likelihood'] = opt_lh  # доп.графа в исх.файле, со ср.знач.вероятностей отображения
            for s, p_lh in lh.items():
                lh_name = str(s) + '_ATOM_PROBABILITY_MAP'
                reaction['meta'][lh_name] = p_lh
                # доп.графы в исх.файле, с вероятностями отображения данного атома реагента на все атомы продукта
            '''
            outputdata.write(reaction)  # запись реакции в исходящий файл
    _, _ = truth(kwargs['input'], kwargs['output'], 0, 0, [], kwargs['debug'])  # Проверка соответствия

