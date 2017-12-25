from multiprocess import Pool, Process, Queue
from pprint import pprint
from timeit import default_timer as timer
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
    print(model['type_model'])
    if model['type_model'] == 'mlp':
        print('hidden_layer_sizes:\t{}\nactivation:\t{}\nsolver:\t{}\nalpha:\t{}'.format(tuple(model['mlp_hls']),
                                                                                         model['mlp_a'],
                                                                                         model['mlp_s'],
                                                                                         model['mlp_alpha']))
    fragger = Fragger(model['min'], model['max'], model['deep'], model['fragment_type'])
    bitstring = Bitstringen(model['bitstring'], model['length'], model['fragment_count'])
    # При предсказании, НИКОГДА не применяем алгоритма Моргана для генерации пар сим./экв. атомов
    pairwise = Pairwise(0, False)

    # dfs_pr = [0, 0, 0]
    total_time = dict.fromkeys(["Bitstring", "Predict", "Munkres", "DFS1", "DFS2", "All"], 0)
    time_ = timer()

    print("Testing set descriptor calculation")
    with open(kwargs['input'], encoding='cp1251') as fr, open(kwargs['output'], 'w') as fw:
        outputdata = RDFwrite(fw)

        for reaction in worker(RDFread(fr), kwargs['debug']):
            p_graph = nx.union_all(reaction['products'])
            reaction['products'] = remap(reaction['products'], {k: k + max(p_graph.nodes()) for k in p_graph.nodes()})

            start_1 = timer()
            time_pr = 0
            y, pairs = [], []
            for x, _, drop_pairs in getXY(reaction, fragger, pairwise, bitstring, model['chunk']):
                '''
                На основании сгенерированного набора битовых строк дескрипторов из обученой модели 
                выводятся значения логорифмов вероятностей проецирования (отображения) атомов
                '''
                pairs.extend(drop_pairs)
                start_2 = timer()
                y.extend(model['model'].predict_log_proba(x))
                end_2 = timer()
                total_time["Predict"] += (end_2 - start_2)
                time_pr += (end_2 - start_2)
            total_time["Bitstring"] += (timer() - start_1 - time_pr)
            _map, map_time = mapping(pairs, y, nx.union_all(reaction['products']), nx.union_all(reaction['substrats']))
            total_time["Munkres"] += map_time[0]
            total_time["DFS1"] += map_time[1]
            total_time["DFS2"] += map_time[2]

            # на основании обученной модели перемаппливаются атомы продукта
            reaction['products'] = remap(reaction['products'], _map)
            outputdata.write(reaction)  # запись реакции в исходящий файл

    total_time["All"] += timer() - time_
    _, _ = truth(kwargs['input'], kwargs['output'], 0, 0, [], kwargs['debug'])  # Проверка соответствия
    if kwargs['debug']:
        pprint(total_time)
