__author__ = 'valiaafo'
import numpy as np

class Prepare(object):
    def __init__(self):
        self.__header_index = 0
        #self.__atomfragcount = {}
        self.__atomfragcount = []
        #self.__atomfragcount_index = 0
        self.__index_bit = []
        self.__bit = []
        #self.__map = {}
        self.__map = []
        #self.__map_index = 0
        self.__y = []
        self.__num_atoms_react = []

    def collect(self, data, header, atomfragcount, mode, reaction):
        if mode == 0:
            atomfragcount = self.__atomfragcount
        else:
            atomfragcount.clear()
        atomfragcount.append({})
        for role,frags in data.items():
            atomfragcount[reaction][role]=[]
            for i,atom in enumerate(frags.keys()):
                atomfragcount[reaction][role].append({})
                atomfragcount[reaction][role][i]['num_a']=atom[0]
                atomfragcount[reaction][role][i]['sym_a']=atom[1]
                atomfragcount[reaction][role][i]['fragments']={}
                for fname,count in frags[atom].items():
                    #print("fname: "+str(fname)+" , count: "+str(count))
                    if mode == 0:
                        if header.get(fname) is None:
                            header[fname] = self.__header_index
                            self.__header_index += 1
                        atomfragcount[reaction][role][i]['fragments'][header[fname]]= count
                    else:
                        if header.get(fname) is not None:
                            atomfragcount[reaction][role][i]['fragments'][header[fname]]= count
            '''
            atomfragcount[reaction][role] = {}
            for atom,fragments in frags.items():
                number_a = atom[0]
                symbol_a = atom[1]
                atomfragcount[reaction][role][number_a] = {}
                atomfragcount[reaction][role][number_a][symbol_a] = {}
                for fname,count in fragments.items():
                    if mode == 0:
                        if header.get(fname) is None:
                            header[fname] = self.__header_index
                            self.__header_index += 1
                        atomfragcount[reaction][role][number_a][symbol_a][header[fname]]= count
                    else:
                        if header.get(fname) is not None:
                            #print("no frag")
                            atomfragcount[reaction][role][number_a][symbol_a][header[fname]]= count
            '''
        #self.__atomfragcount_index +=1
        return atomfragcount,header

    def good_map(self, data, reaction):
        self.__map.append({})
        #print("good_map: ", data)
        for role,molecul in data.items():
            if role != 'meta':
                #chet = []
                #self.__map[reaction][role] = {}
                self.__map[reaction][role] = []
                for m in molecul: #adelia m - данные одной молекуле реагента/продукта
                    for j,atom in enumerate(m['atomlist']):
                        self.__map[reaction][role].append(atom['map'])
                    '''
                    lchet = len(chet)
                    for j,atom in enumerate(m['atomlist']): # j - кол-во атомов в молекуле
                        num_a = j + lchet
                        #self.__map[reaction][role][num_a] = {}
                        chet.append(j)
                        map_a = int(atom['map'])
                        print(map_a)
                        self.__map[role][num_a] = map_a
                    '''
        #self.__map_index += 1
        #maps_dict=self.__map
        #return maps_dict
        maps_list = self.__map
        return maps_list

    def bit_string(self,data,header,mode,type_of_bitstring,j):
        #print("self.__map: ",self.__map)
        index_bit = []
        self.__bit = []
        self.__y = []
        #print("bit_string: ",data)
        for r,react in enumerate(data):#a_deal
            num_r=j+r
            quantity_a=len(react['substrats']) # len(react['substrats'].keys())
            quantity_p=len(react['products'])
            if quantity_a==quantity_p:
                for atom_s in react['substrats']:
                    num_s=atom_s['num_a']
                    sym_s=atom_s['sym_a']
                    frags1=atom_s['fragments']
                    for atom_p in react['products']:
                        sym_p=atom_p['sym_a']
                        if sym_s==sym_p:
                            num_p=atom_p['num_a']
                            frags2=atom_p['fragments']
                            temp_bit_all = []
                            temp_bit_intersection = [0 for x in range(len(header))]
                            com_frag = list(set(frags1).intersection(frags2))
                            for i in range(len(com_frag)):
                                    temp_bit_intersection[com_frag[i]] = 1

                            #common to the three types of
                            if type_of_bitstring == 0:
                                temp_bit_all.extend(temp_bit_intersection)

                            elif type_of_bitstring == 1:
                                temp_bit_s = [0 for x in range(len(header))]
                                for i in range(len(list(frags1.keys()))): #TM changed frags1 to frags
                                    temp_bit_s[list(frags1.keys())[i]] = 1 #TM changed frags1 to frags
                                temp_bit_all.extend(temp_bit_s)
                                temp_bit_p = [0 for x in range(len(header))]
                                for i in range(len(list(frags2.keys()))): #TM changed frags to frags1
                                    temp_bit_p[list(frags2.keys())[i]] = 1 #TM changed frags to frags1
                                temp_bit_all.extend(temp_bit_p)
                                temp_bit_all.extend(temp_bit_intersection)

                            else:
                                temp_bit_s = [0 for x in range(len(header))]
                                a_without_b = list(set(frags1).difference(frags2))
                                for i in range(len(a_without_b)):
                                    temp_bit_s[a_without_b[i]] = 1
                                temp_bit_all.extend(temp_bit_s)
                                temp_bit_p = [0 for x in range(len(header))]
                                b_without_a = list(set(frags2).difference(frags1))
                                for i in range(len(b_without_a)):
                                    temp_bit_p[b_without_a[i]] = 1
                                temp_bit_all.extend(temp_bit_p)
                                temp_bit_all.extend(temp_bit_intersection)

                            #print(temp_bit_all)
                            index_bit.append((num_r,num_s,num_p))
                            self.__bit.append(temp_bit_all)
                            #print(good_map)
                            if mode == 0:
                                if self.__map[num_r]['substrats'][num_s] == self.__map[num_r]['products'][num_p]:
                                    m = 1
                                else:
                                    m = 0
                                self.__y.append(m)
            else:
                print('Error! Number of atoms in reagent and products differ')
        self.__bit = np.array(self.__bit,dtype=bool)
        np.set_printoptions(threshold=np.nan)
        bit_string_res = self.__bit
        if mode == 0:
            y = np.array(self.__y,dtype=bool)
        elif mode == 1:
            y = self.__y
        #print(index_bit)
        return bit_string_res, y, index_bit, header, quantity_a
