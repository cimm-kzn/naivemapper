__author__ = 'stsouko'
import numpy as np


class Prepare(object):
    def __init__(self):
        self.__header = {}
        self.__header_index = 0
        self.__atomfragcount = {}
        self.__atomfragcount_index = 0
        self.__index_bit = []
        self.__bit = []
        self.__map = {}
        self.__map_index = 0
        self.__y = []




    def collect(self, data):
        self.__atomfragcount[self.__atomfragcount_index] = {}
        for role,frags in data.items():
            self.__atomfragcount[self.__atomfragcount_index][role] = {}
            for atom,fragments in frags.items():
                number_a = atom[0]
                symbol_a = atom[1]
                self.__atomfragcount[self.__atomfragcount_index][role][number_a] = {}
                self.__atomfragcount[self.__atomfragcount_index][role][number_a][symbol_a] = {}
                for fname,count in fragments.items():
                    if self.__header.get(fname) is not None:
                        pass
                    else:
                        self.__header[fname] = self.__header_index
                        self.__header_index += 1
                    self.__atomfragcount[self.__atomfragcount_index][role][number_a][symbol_a][self.__header[fname]]= count
        self.__atomfragcount_index +=1
        return self.__atomfragcount

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
                                    # if dat2 == 'element':
                                    #     sym_a = dat_value
                                    #     self.__map[self.__map_index][role][num_a][sym_a] = {}
                                    #     for dat2,dat_value in dat_atom[j].items():
                                            if dat2 == 'map':
                                                map_a = int(dat_value)
                                                self.__map[self.__map_index][role][num_a] = map_a
        self.__map_index += 1
        maps_dict = self.__map
        return maps_dict


    def bit_string(self,data):
        for reaction,react in data.items():
            for role,atom in react.items():
                if role == 'products':
                    pass
                else:
                    for number_a,sym_and_frag in atom.items():
                        num_s = number_a
                        for symbol_a,frags in sym_and_frag.items():
                            sym_s = symbol_a
                            for role,atom in react.items():
                                if role == 'substrats':
                                    pass
                                else:
                                    for number_a,sym_and_frag in atom.items():
                                        num_p = number_a
                                        for symbol_a,frags1 in sym_and_frag.items():
                                            sym_p = symbol_a
                                            if sym_s != sym_p:
                                                # print('no, reaction N',reaction)
                                                # print('atom of substrats N',num_s,', symbol = ',sym_s)
                                                # print('atom of products N',num_p,', symbol = ',sym_p)
                                                # print(fr_s,', ',frags1)
                                                pass
                                            else:
                                                # print('reaction N',reaction)
                                                # print('atom of substrats N',num_s,', symbol = ',sym_s)
                                                # print('atom of products N',num_p,', symbol = ',sym_p)
                                                # print(frags,frags1)
                                                # print(set(frags).intersection(frags1))
                                                temp_bit = [0 for x in range(self.__header_index)]
                                                # print(temp_bit)
                                                com_frag = list(set(frags).intersection(frags1))
                                                # print(com_frag)
                                                for i in range(len(com_frag)):
                                                    temp_bit[com_frag[i]] = 1
                                                    # print(temp_bit)
                                                self.__index_bit.append((reaction,num_s,num_p))
                                                self.__bit.append(temp_bit)
                                                if self.__map[reaction]['substrats'][num_s] == self.__map[reaction]['products'][num_p]:
                                                    m = 1
                                                else:
                                                    m = 0
                                                self.__y.append(m)
        self.__bit = np.array(self.__bit)
        np.set_printoptions(threshold=np.nan)
        # print(self.__bit,type(self.__bit),type(self.__bit[0]),type(self.__bit[0][0]))
        bit_string_res = self.__bit
        y = np.array(self.__y)
        index_bit = self.__index_bit
        return bit_string_res, y,index_bit









