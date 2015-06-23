#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of naivemapper.
#
# naivemapper is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import argparse

from core.fragger import Fragger
from core.RDFread import RDFread
from core.version import version
from core.RDFwrite import RDFwrite
from core.Prepare import Prepare
from core.Models import Models

import pickle

import itertools

from sklearn.naive_bayes import BernoulliNB


def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield dict(chunk)

def main():
    rawopts = argparse.ArgumentParser(description="Naive Reactions Automapper",
                                      epilog="Copyright 2015 Ramil Nugmanov <stsouko@live.ru>")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default='input.rdf', help="input RDF ")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="output RDF")
    rawopts.add_argument("--min", "-m", type=int, default=1, help="minimal fragments length")
    rawopts.add_argument("--max", "-M", type=int, default=8, help="maximal fragments length")
    rawopts.add_argument("--model", "-n", type=str, default="model", help="name of file with model")
    rawopts.add_argument("--predict", "-p", type=bool, default=1, help="mode of the program: 0 - learning, 1 - prediction") # default=0 !!!
    rawopts.add_argument("--bitstring", "-b", type=int, default=2, help="type of united bitstring for atomic bitstrings A and B: 0 - A*B, 1- A+B+A*B, 2 - A!B + B!A + A*B")
    options = vars(rawopts.parse_args())

    inp = RDFread(options['input'])
    if not inp.chkRDF():
        print('rdf incorrect')
        return 0

    bit = []
    atomfragcount = {}
    bitstring = options['bitstring']

    if options['predict'] == 0:
        fragger = Fragger(**options)
        collector = Prepare()
        e = 0
        header = {}
        Model = Models()
        model = BernoulliNB()

        for i, data in enumerate(inp.readdata()):
                        if i % 100 == 0 and i:
                            print("reaction: %d" % (i + 1))
                        #res = calc.firstcgr(data)
                        try:
                            res = fragger.get(data)
                            maps_dict = collector.good_map(data)
                            new_data,header = collector.collect(res,header,atomfragcount,0)

                        except:
                            e += 1
                            print("Error: %d" % (i + 1))

        for chunk in grouper(10,new_data.items()):
            print('chunk = ',chunk)
            bit,y,index_bit,header,quantity_a = collector.bit_string(chunk,header,0,bitstring)
            print('Modeling in progress... X size: ', len(index_bit))
            model = Model.learning(bit,y,model)
        print(y, index_bit)
        filename = options['model']
        folder = 'trained_models\\'
        filename_extension = '.pickle'
        filename = folder+filename+filename_extension


        model_and_header = {'model':model,'header':header}

        with open(filename,'wb') as f_tr_model_and_header:
            pickle.dump(model_and_header,f_tr_model_and_header)
        f_tr_model_and_header.close()


    else:
        out = RDFwrite(options['output'])
        fragger = Fragger(**options)
        collector = Prepare()
        e = 0

        filename = options['model']
        folder = 'trained_models\\'
        filename_extension = '.pickle'
        filename = folder+filename+filename_extension

        Model = Models()

        with open(filename,'rb') as f_tr_model_and_header:
            model_and_header = pickle.load(f_tr_model_and_header)
        f_tr_model_and_header.close()
        header = model_and_header['header']
        model = model_and_header['model']

        for i, data in enumerate(inp.readdata()):
            if i % 100 == 0 and i:
                print("reaction: %d" % (i + 1))
            #res = calc.firstcgr(data)

            # try:
            res = fragger.get(data)
            new_data,header = collector.collect(res,header,atomfragcount,1)
            bit,y,index_bit,header,quantity_a = collector.bit_string(new_data,header,1,bitstring)
            probabilities = Model.predict(model,bit)
            index_map = Model.mapping(index_bit,probabilities,quantity_a)
            ind = 1
            s_map = [(x,z) for x,y in enumerate(data['substrats']) for z in range(len(y['atomlist']))] # порядковый номер соответствует номеру атома в index_map, в кортеже первый номер соответствует номеру молекулы, второй - номеру атома в молекуле
            p_map = [(x,z) for x,y in enumerate(data['products']) for z in range(len(y['atomlist']))]

            for role,dat in data.items():
                if role == 'substrats':
                    for i,datum in enumerate(dat):
                        for list_type, dat_atom in datum.items():
                            if list_type == 'atomlist':
                                for j,properties in enumerate(dat_atom):
                                    for prop,count in properties.items():
                                        if prop == 'map':
                                            data['substrats'][i]['atomlist'][j]['map'] = ind
                                            l = s_map.index((i,j)) # получили положение кортежа в списке реагентов, соответствующее сквозному номеру атома реагента
                                            k = index_map[l] # по сквозному номеру атома реагента находим атом продукта, которому он соответствует
                                            m,n = p_map[k] # получаем номер молекулы и номер атома продукта в данной молекуле
                                            data['products'][m]['atomlist'][n]['map'] = ind
                                            ind += 1
            out.write(data)

            # except:
            #     e += 1
            #     print("Error: %d" % (i + 1))
            #out.write(data)

    # print(model.class_count_, model.class_log_prior_,model.feature_count_[0],model.feature_count_[1])
    print("Checked %d reactions. %d reactions consist exception errors" % (i + 1, e))
    return 0


if __name__ == '__main__':
    main()
