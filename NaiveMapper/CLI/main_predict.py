from pprint import pprint
from random import shuffle
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
        tmp.append(graph.remap({k: v for k, v in maps.items() if k in graph}, copy=True))
    return tmp


def predict_core(**kwargs):
    # Загружаем уже существующую модель в Наивный Бейсовский классификатор
    model = pickle.load(kwargs['model'])
    fragger = Fragger(f_type=model['fragment_type'], _min=model['min'], _max=model['max'], _deep=model['deep'],
                      _min2=model['min2'], _max2=model['max2'], _fuzzy=model['fuzzy'])
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
            p_nodes = [k + max(p_graph.nodes()) for k in p_graph.nodes()]
            shuffle(p_nodes)
            reaction.products._MindfulList__data = remap(reaction['products'], {k: v for k, v in zip(p_graph.nodes(), p_nodes)})

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
            reaction.products._MindfulList__data = remap(reaction['products'], _map)
            outputdata.write(reaction)  # запись реакции в исходящий файл

    total_time["All"] += timer() - time_
    _, _ = truth(kwargs['input'], kwargs['output'], 0, 0, [], kwargs['debug'])  # Проверка соответствия
    if kwargs['debug']:
        pprint(total_time)
