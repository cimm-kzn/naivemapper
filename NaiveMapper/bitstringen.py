#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of naivemapper.
#
#  naivemapper is free software; you can redistribute it and/or modify
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
import pandas as pd
import hashlib
from itertools import product
from typing import Dict, Tuple, Iterable, List


class Bitstringen(object):
    def __init__(self, b_type, b_length, f_count=False):
        self.__combiner = self.combo1 if b_type == 1 else self.combo2 if b_type == 2 else self.combo3 if b_type == 3 \
            else self.combo4 if b_type == 4 else self.combo5 if b_type == 5 else self.combo0
        self.__b_type = b_type

        self.__length = b_length
        self.__type = f_count

    def get(self, substrats: Dict[int, Dict[str, int]], products: Dict[int, Dict[str, int]],
            pairs: Iterable[Tuple[int, int]]) -> pd.DataFrame:
        """
        из словарей реагентов и продуктов, где ключами являются номер атома, а значением словарь типа:
        {картеж_фрагмента: кол-во}, для каждого атома реагента/продукта берется словарь фрагментов(дискрипторы),
        который обрабытываются функцией get_bitstring, для каждой пары атом_реагент и атом_продукт комбинируются,
        в зависимости от b_type, желаемая битовая строка.
        На выходе имеется таблица, элементами(строками) которой являются битовые строки для каждой пары атомов.
        """
        s_cache, p_cache = {}, {}
        combos = []
        for s_atom, p_atom in pairs:
            if self.__type:  # сли необходимо учитывать кол-во фрагментов каждого типа
                # Перечисление фрагментов конкатенируется с его названием
                s_set = s_cache[s_atom] if s_atom in s_cache \
                    else s_cache.setdefault(s_atom, set([str(j) + k for k, v in substrats[s_atom].items()
                                                         for j in range(1, v + 1)]))
                p_set = p_cache[p_atom] if p_atom in p_cache \
                    else p_cache.setdefault(p_atom, set([str(j1) + k1 for k1, v1 in products[p_atom].items()
                                                         for j1 in range(1, v1 + 1)]))
            else:  # Игнорируется информация о количестве фрагментов у атомов
                s_set = s_cache[s_atom] if s_atom in s_cache else s_cache.setdefault(s_atom, substrats[s_atom].keys())
                p_set = p_cache[p_atom] if p_atom in p_cache else p_cache.setdefault(p_atom, products[p_atom].keys())
            combos.append(self.__combiner(s_set, p_set))

        return pd.DataFrame(combos)

    def get_bitstring(self, descriptors: List) -> pd.Series:
        """
        Дискриптор(словарь вида {(фрагмент):кол-во}).
        Фрагменты атома реагента/продукта ХЕШируется(т.е. преобритается уникальный код).
        Элемент битовой строки имеет значение 1 если такой фрагмент от данного атома существует, и 0 - в ином случае.
        """
        bitstring = pd.Series(0, index=range(self.__length))
        for k in descriptors:
            hs = int(hashlib.md5(k.encode()).hexdigest(), 16)
            bitstring[hs % self.__length] = 1
            bitstring[hs // self.__length % self.__length] = 1
            # l = bitstring[bitstring.isin([1])].index.tolist()
        return bitstring > 0

    def combo0(self, s_set, p_set):
        # A*B (в битовой строке ставится 1 если фрагменты присутствовали в реагенте И в продукте)
        list_frag = list(x for x in s_set if x in p_set)
        return self.get_bitstring(list_frag)

    def combo1(self, s_set, p_set):
        """
        A + B + A*B (битовая строка утраивается.
        Состоит из бит.строк реагента, продукта и их умножения(в продукте И в реагенте)
        """
        return pd.concat([self.get_bitstring(s_set), self.get_bitstring(p_set), self.combo0(s_set, p_set)],
                         keys=['A', 'B', 'A*B'])

    def combo2(self, s_set, p_set):
        """
        A!B + B!A + A*B (утроенная битовая строка), состоящая из бит.строк:
        - есть в реагенте НО НЕ в продукте,
        - есть в продукте НО НЕ в реагенте,
        - есть в реагенте И в продукте.
        """
        only_a = set(x for x in s_set if x not in p_set)
        only_b = set(y for y in p_set if y not in s_set)
        return pd.concat([self.get_bitstring(only_a), self.get_bitstring(only_b), self.combo0(s_set, p_set)],
                         keys=['A!B', 'B!A', 'A*B'])

    def combo3(self, s_set, p_set):
        """
        A + B (удвоенная битовая строка), состоящая из бит.строк:
        - есть в реагенте,
        - есть в продукте.
        Работает ОЧЕНЬ плохо!
        """
        return pd.concat([self.get_bitstring(s_set), self.get_bitstring(p_set)], keys=['A', 'B'])

    def combo4(self, s_set, p_set):
        """
        (A xor B) + A*B (удвоенная битовая строка), состоящая из бит.строк:
        - если бит присутствует в реагенте или в продукте, но не в обоих сразу,
        - есть в реагенте И в продукте."""
        return pd.concat([self.get_bitstring(s_set ^ p_set), self.combo0(s_set, p_set)], keys=['A xor B', 'A*B'])

    def combo5(self, s_set, p_set):
        """
        Паросочетания всех фрагментов реагентов и продуктов
        """
        return self.get_bitstring([s_fr + '&' + p_fr for s_fr, p_fr in product(s_set, p_set)])
