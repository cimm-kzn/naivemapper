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
import pickle
import itertools

from core.fragger import Fragger
from core.RDFread import RDFread1
from core.version import version
from core.RDFwrite import RDFwrite
from core.Prepare import Prepare
from core.Models import Models
from core.Properties import Properties

from CGRtools.CGRcore import CGRcore
from CGRtools.FEAR import FEAR
from CGRtools.RDFread import RDFread
from CGRtools.SDFwrite import SDFwrite

from sklearn.naive_bayes import BernoulliNB

def grouper(n, iterable):
    it = iter(x for x in iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           break
       yield list(chunk)

def main():
    rawopts = argparse.ArgumentParser(description="Naive Bayes Reactions Automapper",
                                      epilog="Copyright 2015 Ramil Nugmanov <stsouko@live.ru>")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default='input.rdf', help="input RDF ")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="output RDF")
    rawopts.add_argument("--min", "-m", type=int, default=1, help="minimal fragments length")
    rawopts.add_argument("--max", "-M", type=int, default=8, help="maximal fragments length")
    rawopts.add_argument("--model", "-n", type=str, default="model", help="name of file with model")
    rawopts.add_argument("--predict", "-p", type=int, default=0,
                         help="mode of the program: 0 - learning, 1 - prediction, 2 - cross_validation") # default=0 !!!
    rawopts.add_argument("--bitstring", "-b", type=int, default=2,
                         help="type of united bitstring for atomic bitstrings A and B: 0-A*B, 1-A+B+A*B, 2-A!B+B!A+A*B")
    options = vars(rawopts.parse_args())
    inp = RDFread1(options['input'])
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
                #print(res)
                #res = role: (atomNo, atomSymb, atomCharge?): frag: fragCount
                maps_dict = collector.good_map(data) #похоже не нужна
                new_data,header = collector.collect(res,header,atomfragcount,0,i)
                # newdata is dictionary with atoms and its fragments, reacNo: role: atomNo: atomSymbol: frag: fragCount
                # header is dictionary of fragmentDescription (tuple): fragmentNo
            except:
                e += 1
                print("Error: %d" % (i + 1))
        #print("newdata", new_data)
        print("Totally number of data for train:", len(new_data))
        model, header = Properties.train(new_data, collector, bitstring, Model, model, header)
        model_and_header = {'model':model,'header':header}
        #print(y, index_bit)

        filename = options['model']
        folder = 'trained_models/'
        filename_extension = '.pickle'
        filename = folder+filename+filename_extension
        print(filename)
        with open(filename,'wb') as f_tr_model_and_header:
            pickle.dump(model_and_header,f_tr_model_and_header)
        f_tr_model_and_header.close()


    elif options['predict'] == 1:
        out = RDFwrite(options['output'])
        fragger = Fragger(**options)
        collector = Prepare()
        e = 0

        filename = options['model']
        folder = 'trained_models/'
        filename_extension = '.pickle'
        filename = folder+filename+filename_extension

        Model = Models()

        with open(filename,'rb') as f_tr_model_and_header:
            model_and_header = pickle.load(f_tr_model_and_header)
        f_tr_model_and_header.close()
        header = model_and_header['header']
        model = model_and_header['model']
        final_data = {}
        for i, data in enumerate(inp.readdata()):
            if i % 10 == 0:
                print("reaction: %d" % (i + 1))
            #res = calc.firstcgr(data)

            # try:
            res = fragger.get(data)
            new_data,header = collector.collect(res,header,atomfragcount,1, i)
            bit,y,index_bit,header,quantity_a = collector.bit_string(new_data,header,1,bitstring)
            probabilities = Model.predict(model,bit)
            index_map = Model.mapping(index_bit,probabilities,quantity_a)
            s_map = [(x,z) for x,y in enumerate(data['substrats']) for z in range(len(y['atomlist']))]
            ''' порядковый номер соответствует номеру атома в index_map,
            в кортеже первый номер соответствует номеру молекулы, второй - номеру атома в молекуле'''
            p_map = [(x,z) for x,y in enumerate(data['products']) for z in range(len(y['atomlist']))]
            #print("data", data)
            data = Properties.test_write(data, s_map, p_map, index_map, options)
            out.write(data)
            final_data[i] = data
        #out.write(final_data)
        print("final_data", final_data)

        #out.write(final_data)

    elif options['predict'] == 2:
        e = 0
        #header = {}
        N_folds = 5
        #data_all={}
        data_all = []

        CGR=CGRcore(type='0', stereo=False, b_templates=None, balance=0, c_rules= None, e_rules=None)
        fear=FEAR()
        #SDF=SDFwrite()

        # parse files

        for i, data in enumerate(inp.readdata()):
            if i % 100 == 0 and i:
                print("reaction: %d" % (i + 1))
            data_all.append(data)

        #generates folds

        correct_mapping_CV = 0
        correct_mapping_Fit = 0
        incorrect_mapping_CV = 0
        incorrect_mapping_Fit = 0
        correct_CV=[]
        correct_Fit=[]
        incorrect_CV = []
        incorrect_Fit=[]

        #data_predicted = {}
        for fold in range(N_folds+1):
            bit = []
            #atomfragcount = {}
            atomfragcount = []
            bitstring = options['bitstring']
            header = {}
            Model = Models()
            model = BernoulliNB()
            fragger = Fragger(**options)
            collector = Prepare()

            #data_train = {}
            #data_test = {}
            data_train = []
            data_test = []

            if fold != N_folds: # for cross-validation performance estimation
                print("Fold", fold+1, "/", N_folds)
                for i in range(len(data_all)):
                    if i % N_folds != fold:
                        data_train.append(data_all[i])
                    else:
                        data_test.append(data_all[i])
            else: # for fitting performance estimation and final model preparation
                print("All data")
                for i in range(len(data_all)):
                    data_train.append(data_all[i])
                    data_test.append(data_all[i])

            #training stage

            print("Training set descriptor calculation")
            #for j, data in data_train.items():
            for j,data in enumerate(data_train):
                #print("j: "+str(j)+", data: "+ str(data))
                try:
                    res = fragger.get(data)
                    maps_list=collector.good_map(data,j)
                    new_data, header=collector.collect(res, header, atomfragcount, 0, j)
                    # res is dictionary role: atom(nomber,'symbol','zaryad'): frag: fragCount
                    # maps_dict is dictionary test-set reacNo: role: atomNo: (map)
                    # (old)newdata is dictionary with atoms and its fragments, reacNo: role: atomNo: atomSymbol: frag: fragCount
                    # (new)newdata is list [j(=reacNo) {role [{'num_a', 'sym_a', 'fragments': frag: fragCount},...]}]
                    # header is dictionary of fragmentDescription (tuple): fragmentNo
                except:
                    e += 1
                    print("Error: %d" % (j + 1))

            print(len(data_train), "reactions ready")

            model, header = Properties.train(new_data, grouper, collector, bitstring, Model, model, header)
            model_and_header = {'model':model,'header':header}

            filename = options['model']
            folder = 'trained_models/'
            filename_extension = '.pickle'
            filename = folder+filename+'_fold'+str(fold)+filename_extension
            print("filename:", filename)
            with open(filename,'wb') as f_tr_model_and_header:
                pickle.dump(model_and_header,f_tr_model_and_header)
            f_tr_model_and_header.close()

            # prediction stage

            print("Testing set descriptor calculation")
            #data_predicted = {}
            data_predicted = []
            filename_pred = "output_fold"+str(fold)+"pred.rdf"
            filename_test = "output_fold"+str(fold)+"test.rdf"
            print("filename_pred: ", filename_pred, "filename_test: ", filename_test)
            out_pred = RDFwrite(filename_pred)
            out_test = RDFwrite(filename_test)

            #for j, data in data_test.items():
            for j,data in enumerate(data_test):
                try:
                    res = fragger.get(data)
                    new_data2,header = collector.collect(res,header,atomfragcount,1,j)
                    bit,y,index_bit,header,quantity_a = collector.bit_string(new_data2,maps_list,header,1,bitstring,j)
                except:
                    e += 1
                    print("Error: %d" % (j + 1))

                probabilities = Model.predict(model,bit)
                index_map = Model.mapping(index_bit,probabilities,quantity_a)
                s_map = [(x,z) for x,y in enumerate(data['substrats']) for z in range(len(y['atomlist']))]
                ''' порядковый номер соответствует номеру атома в index_map, в кортеже первый номер
                соответствует номеру молекулы, второй - номеру атома в молекуле'''
                p_map = [(x,z) for x,y in enumerate(data['products']) for z in range(len(y['atomlist']))]
                #print("data", data)
                out_test.write(data)
                data = Properties.test_write(data, s_map, p_map, index_map, options)
                out_pred.write(data)
                data_predicted.append(data)
                #print("j "+str(j)+", id "+data['meta']['CdId'])
            out_pred.close()
            out_test.close()
            #print(len(data_predicted))

            # comparison stage - valid for balanced reactions only

            filename_SDF = "out_fold"+str(fold)+".sdf"
            SDF = SDFwrite(open(filename_SDF, 'w'))

            with open(filename_pred) as predfile, open(filename_test) as testfile:
                reac=1
                for pred,test in zip(RDFread(predfile).readdata(),RDFread(testfile).readdata()):
                    predCGR=CGR.getCGR(pred)
                    testCGR=CGR.getCGR(test)
                    predHash=fear.getreactionhash(predCGR)
                    testHash=fear.getreactionhash(testCGR)

                    if fold != N_folds:
                        if predHash==testHash:
                            correct_mapping_CV+=1
                            correct_CV.append(reac)
                        else:
                            incorrect_mapping_CV+=1
                            incorrect_CV.append(reac)
                        reac+=1
                    else:
                        if predHash==testHash:
                            correct_mapping_Fit+=1
                            correct_Fit.append(reac)
                        else:
                            incorrect_mapping_Fit+=1
                            incorrect_Fit.append(reac)
                        reac+=1
                    #reac+=1
                    tmp=CGR.getformattedcgr(predCGR)
                    tmp['meta']={}
                    SDF.writedata(tmp)
                print('corr_CV:', correct_CV, 'incorr_CV:', incorrect_CV)
                correct_CV.clear()
                incorrect_CV.clear()
            print("percent correct_CV: ", str(100*correct_mapping_CV/(correct_mapping_CV+incorrect_mapping_CV))+"%")
        print('corr_Fit:', correct_Fit, 'incorr_Fit:', incorrect_Fit)
        print("percent correct_Fit: ",  str(100*correct_mapping_Fit/(correct_mapping_Fit+incorrect_mapping_Fit))+"%")

    # print(model.class_count_, model.class_log_prior_,model.feature_count_[0],model.feature_count_[1])
    print("Checked %d reactions. %d reactions consist exception errors" % (i + 1, e))
    return 0

if __name__ == '__main__':
    main()
