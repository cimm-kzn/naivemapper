__author__ = 'valiaafo'
import numpy as np


class Prepare(object):
    def __init__(self):
        #self.__header = {}
        self.__header_index = 0
        self.__atomfragcount = {}
        self.__atomfragcount_index = 0
        self.__index_bit = []
        self.__bit = []
        self.__map = {}
        self.__map_index = 0
        self.__y = []
        self.__num_atoms_react = []




    def collect(self, data, header,atomfragcount,mode):
        if mode == 0:
            atomfragcount = self.__atomfragcount
            # print('atomfragcount from collect = ',atomfragcount)
        else:
            atomfragcount.clear()
        atomfragcount[self.__atomfragcount_index] = {}
        for role,frags in data.items():
            atomfragcount[self.__atomfragcount_index][role] = {}
            for atom,fragments in frags.items():
                number_a = atom[0]
                symbol_a = atom[1]
                atomfragcount[self.__atomfragcount_index][role][number_a] = {}
                atomfragcount[self.__atomfragcount_index][role][number_a][symbol_a] = {}
                for fname,count in fragments.items():
                    if mode == 0:
                        if header.get(fname) is not None:
                            pass
                        else:
                            header[fname] = self.__header_index
                            self.__header_index += 1
                        atomfragcount[self.__atomfragcount_index][role][number_a][symbol_a][header[fname]]= count
                    else:
                        # if header.get(fname) is None:
                        #     pass
                        # else:
                            atomfragcount[self.__atomfragcount_index][role][number_a][symbol_a][header[fname]]= count

        self.__atomfragcount_index +=1
        # print('atomfragcount from collect 2 = ',atomfragcount)
        return atomfragcount,header

    def good_map(self, data):
        self.__map[self.__map_index] = {}
        for role,dat in data.items():
            if role == 'meta':
                pass
            else:
                chet = []
                self.__map[self.__map_index][role] = {}
                for i in range(len(dat)):
                    for list_type, dat_atom in dat[i].items():
                        if list_type == 'atomlist':
                            lchet = len(chet)
                            for j in range(len(dat_atom)):
                                num_a = j + lchet
                                self.__map[self.__map_index][role][num_a] = {}
                                chet.append(j)
                                for dat2,dat_value in dat_atom[j].items():
                                    if dat2 == 'map':
                                       map_a = int(dat_value)
                                       self.__map[self.__map_index][role][num_a] = map_a
        self.__map_index += 1
        maps_dict = self.__map
        # print(maps_dict)
        return maps_dict


    def bit_string(self,data,header,mode):
        print('data = ',data)
        index_bit = []
        if mode == 1:
            self.__bit = []
            self.__bit.clear()
        quantity_a = 0
        for reaction,react in data.items():
            for role,atom in react.items():
                if role == 'products':
                    pass
                else:
                    for number_a,sym_and_frag in atom.items():
                        num_s = number_a
                        if num_s > quantity_a:
                            quantity_a +=1
                        for symbol_a,frags in sym_and_frag.items():
                            # print('frags = ',frags)
                            sym_s = symbol_a
                            for role,atom in react.items():
                                if role == 'substrats':
                                    pass
                                else:
                                    for number_a,sym_and_frag in atom.items():
                                        num_p = number_a
                                        for symbol_a,frags1 in sym_and_frag.items():
                                            # print('frags1 = ',frags1)
                                            sym_p = symbol_a
                                            if sym_s != sym_p:
                                                # print('no, reaction N',reaction)
                                                # print('atom of substrats N',num_s,', symbol = ',sym_s)
                                                # print('atom of products N',num_p,', symbol = ',sym_p)
                                                pass
                                            else:
                                                # print('reaction N',reaction)
                                                # print('atom of substrats N',num_s,', symbol = ',sym_s)
                                                # print('atom of products N',num_p,', symbol = ',sym_p)
                                                # print(frags,frags1)
                                                # print(set(frags).intersection(frags1))
                                                temp_bit = [0 for x in range(len(header))]
                                                # print(temp_bit)
                                                com_frag = list(set(frags).intersection(frags1))
                                                # print('com_frag = ',com_frag)
                                                for i in range(len(com_frag)):
                                                    temp_bit[com_frag[i]] = 1
                                                    # print('temp_bit = ',temp_bit)
                                                index_bit.append((reaction,num_s,num_p))
                                                # print('self.__index_bit = ',self.__index_bit)
                                                # print('self.__bit 1 = ',self.__bit)
                                                self.__bit.append(temp_bit)
                                                # print('self.__bit = ',self.__bit)
                                                if mode == 0:
                                                    if self.__map[reaction]['substrats'][num_s] == self.__map[reaction]['products'][num_p]:
                                                        m = 1
                                                    else:
                                                        m = 0
                                                    self.__y.append(m)
        self.__bit = np.array(self.__bit)
        # print('np.array(self.__bit) = ',self.__bit)
        np.set_printoptions(threshold=np.nan)
        bit_string_res = self.__bit
        # print('bit_string_res = ',bit_string_res)
        if mode == 0:
            y = np.array(self.__y)
        else:
            y = self.__y
        # index_bit = self.__index_bit
        # print('index_bit = ',index_bit)
        # print('y = ',y)
        # print('header = ',header)
        # print(len(bit_string_res))
        quantity_a += 1
        # print('kol = ',quantity_a)
        return bit_string_res,y,index_bit,header,quantity_a









