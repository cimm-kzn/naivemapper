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
    options = vars(rawopts.parse_args())

    inp = RDFread(options['input'])
    if not inp.chkRDF():
        print('rdf incorrect')
        return 0

    bit = []
    atomfragcount = {}


    if options['predict'] == 0:
        fragger = Fragger(**options)
        collector = Prepare()
        e = 0
        header = {}


        for i, data in enumerate(inp.readdata()):
            print(data)
            if i % 100 == 0 and i:
                print("reaction: %d" % (i + 1))
            #res = calc.firstcgr(data)
            try:
                res = fragger.get(data)
                print('res = ',res)
                maps_dict = collector.good_map(data)
                print(maps_dict)
                new_data,header = collector.collect(res,header,atomfragcount,0)
                print(header)

            except:
                e += 1
                print("Error: %d" % (i + 1))

        bit,y,index_bit,header,quantity_a = collector.bit_string(new_data,header,0)
        print(bit)
        print(len(bit))
        print(y)

        filename = options['model']
        folder = 'trained_models\\'
        filename_extension = '.pickle'
        filename = folder+filename+filename_extension

        Model = Models()
        new_model = Model.learning(bit,y)
        model_and_header = {'model':new_model,'header':header}

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
            print(data)
            if i % 100 == 0 and i:
                print("reaction: %d" % (i + 1))
            #res = calc.firstcgr(data)
            try:
                res = fragger.get(data)
                print('res = ',res)
                new_data,header = collector.collect(res,header,atomfragcount,1)
                print('new_data = ',new_data)
                bit,y,index_bit,header,quantity_a = collector.bit_string(new_data,header,1) #,bit
                print('bit',bit)
                print('type of bit',type(bit))
                print('index_bit = ',index_bit)
                probabilities = Model.predict(model,bit)
                print('probabilities = ',probabilities)
                print('len of probabilities =',len(probabilities))
                index_map = Model.mapping(index_bit,probabilities,quantity_a)
                ind = 1
                for role,dat in data.items():
                    if role == 'substract':
                        for i in range(len(dat)):
                            for list_type, dat_atom in dat[i].items():
                                if list_type == 'atomlist':
                                    for j in range(len(dat_atom)):
                                        for prop,count in dat_atom[j].items():
                                            if prop == 'map':
                                                print(data[role][i][list_type][j][prop])
                                                data[role][i][list_type][j][prop] = ind
                                                ind += 1
                    elif role == 'product':
                        for i in range(len(dat)):
                            for list_type, dat_atom in dat[i].items():
                                if list_type == 'atomlist':
                                    for j in range(len(dat_atom)):
                                        for prop,count in dat_atom[j].items():
                                            if prop == 'map':
                                                pass



            except:
                e += 1
                print("Error: %d" % (i + 1))



    print("Checked %d reactions. %d reactions consist exception errors" % (i + 1, e))
    return 0


if __name__ == '__main__':
    main()
