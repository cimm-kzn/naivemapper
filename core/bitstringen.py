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
from typing import Dict, Tuple, Iterable


class Bitstringen(object):
    def __init__(self, b_type=0, b_length=1024):
        self.__combiner = self.combo1 if b_type == 1 else self.combo2 if b_type == 2 else self.combo0
        self.__length = b_length

    def get(self, substrats: Dict[int, Dict[Tuple, int]], products: Dict[int, Dict[Tuple, int]],
            pairs: Iterable[Tuple[int, int]]) -> pd.DataFrame:
        """
        из словарей реагентов и продуктов, где ключами являются номер атома, а значением словарь(картеж_фрагмента: кол-во),
        для каждого атома реагента/продукта берется словарь фрагментов(дискрипторы), который обрабытываются функцией get_bitstring,
        для каждой пары атом_реагент и атом_продукт комбинируются, в зависимости от b_type, желаемая битовая строка.
        На выходе имеется таблица, элементами(строками) которой являются битовые строки для каждой пары атомов.
        """
        s_cache, p_cache = {}, {}
        combos = []
        for s_atom, p_atom in pairs:
            s_bitstring = s_cache[s_atom] if s_atom in s_cache else s_cache.setdefault(s_atom, self.get_bitstring(substrats[s_atom]))
            p_bitstring = p_cache[p_atom] if p_atom in p_cache else p_cache.setdefault(p_atom, self.get_bitstring(products[p_atom]))
            combos.append(self.__combiner(s_bitstring, p_bitstring))
        return pd.DataFrame(combos)

    def get_bitstring(self, descriptors: Dict[Tuple, int]) -> pd.Series:
        """
        Дискриптор(словарь вида {(фрагмент):кол-во}). Фрагменты атома реагента/продукта ХЕШируется(т.е. преобритается уникальный код)
        Элемент битовой строки имеет значение 0 если такого фрагмента от данного атома не существует, и 1 - в ином случае.
        """
        bitstring = pd.Series(0, index=range(self.__length))
        for k in descriptors:
            hs = hash(k)
            bitstring[hs % self.__length] = 1
            bitstring[hs // self.__length % self.__length] = 1
        return bitstring

    def combo0(self, s_atom: pd.Series, p_atom: pd.Series) -> pd.Series:
        # A*B (в битовой строке ставится 1 если фрагменты присутствовали в реагенте И в продукте)
        return (s_atom*p_atom) > 0

    def combo1(self, s_atom: pd.Series, p_atom: pd.Series) -> pd.Series:
        # A + B + A*B (битовая строка утраивается. Состоит из бит.строк реагента, продукта и их умножения(в продукте И в реагенте)
        return pd.concat([s_atom, p_atom, s_atom*p_atom], keys=['A', 'B', 'A*B']) > 0

    def combo2(self, s_atom: pd.Series, p_atom: pd.Series) -> pd.Series:
        # A!B + B!A + A*B (утроенная битовая строка, состоящая из бит.строк:
        # есть в реагенте НО НЕ в продукте, есть в продукте НО НЕ в реагенте, и есть в реагенте И в продукте
        return pd.concat([s_atom-p_atom, p_atom-s_atom, s_atom*p_atom], keys=['A-B', 'B-A', 'A*B']) > 0
