#!/usr/bin/env python
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
from collections import defaultdict
import os
import subprocess as sp
import numpy


def main():
    print("This file is part of naivemapper.")
    return 0


class RDFread(object):
    def __init__(self, file):
        self.__RDFfile = file

    def chkRDF(self):
        try:
            f = open(self.__RDFfile)
            if "$RDFILE 1" in f.next():
                return True
            else:
                return False
        except IOError:
            return False

    def readdata(self):
        """
        парсер RDF файлов.
        """
        with open(self.__RDFfile) as f:
            ir = -1
            im = -1
            atomcount = -1
            bondcount = -1
            failkey = True
            reaction = None
            meta = None
            rid = 0
            for n, line in enumerate(f):
                if not failkey and "$RXN" not in line[0:4]:
                    print rid
                    continue
                elif "$RXN" in line[0:4]:
                    rid += 1
                    if reaction:
                        yield reaction
                    reaction = {'substrats': [], 'products': [], 'meta': defaultdict(str)}
                    meta = None
                    ir = n + 4
                    failkey = True
                elif n == ir:
                    try:
                        substrats, products = int(line[0:3]), int(line[3:6])
                    except:
                        failkey = False
                        reaction = None
                elif "$MOL" in line[0:4]:
                    molecule = {'atomlist': [], 'bondmatrix': None, 'bondlist': []}
                    im = n + 4
                elif n == im:
                    try:
                        atomcount = int(line[0:3]) + im
                        bondcount = int(line[3:6]) + atomcount
                        if atomcount == bondcount:
                            molecule['bondmatrix'] = numpy.zeros((1, 1), dtype=float)
                    except:
                        failkey = False
                        reaction = None
                elif n <= atomcount:
                    molecule['atomlist'].append(dict(element=line[31:34].strip(), izotop=line[34:36].strip(),
                                                     charge=int(line[38:39]), map=int(line[60:63]),
                                                     x=float(line[0:10]), y=float(line[10:20]), z=float(line[20:30])))
                elif n <= bondcount:
                    try:
                        if molecule['bondmatrix'] is None:
                            mlen = len(molecule['atomlist'])
                            molecule['bondmatrix'] = numpy.zeros((mlen, mlen), dtype=float)
                        a1, a2 = int(line[0:3]) - 1, int(line[3:6]) - 1
                        molecule['bondlist'].append((a1, a2, int(line[6:9]), int(line[9:12])))
                    except:
                        failkey = False
                        reaction = None
                elif "M  END" in line:
                    try:
                        if len(reaction['substrats']) < substrats:
                            reaction['substrats'].append(molecule)
                        else:
                            reaction['products'].append(molecule)
                    except:
                        failkey = False
                        reaction = None
                elif '$DTYPE' in line:
                    meta = line[7:].strip()
                elif '$RFMT' not in line and meta:
                    reaction['meta'][meta] += line.strip("$DATUM").strip() + ' '
            else:
                if reaction:
                    yield reaction


if __name__ == '__main__':
    main()
