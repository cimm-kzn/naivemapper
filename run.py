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


def main():
    rawopts = argparse.ArgumentParser(description="Naive Reactions Automapper",
                                      epilog="Copyright 2015 Ramil Nugmanov <stsouko@live.ru>")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default='input.rdf', help="input RDF ")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="output RDF")
    rawopts.add_argument("--min", "-m", type=int, default=1, help="minimal fragments length")
    rawopts.add_argument("--max", "-M", type=int, default=8, help="maximal fragments length")
    options = vars(rawopts.parse_args())

    inp = RDFread(options['input'])
    if not inp.chkRDF():
        print('rdf incorrect')
        return 0

    out = RDFwrite(options['output'])
    fragger = Fragger(**options)
    collector = Prepare()
    e = 0

    for i, data in enumerate(inp.readdata()):
        if i % 100 == 0 and i:
            print("reaction: %d" % (i + 1))
        #res = calc.firstcgr(data)
        try:
            res = fragger.get(data)
            print(res)
            collector.collect(res)
        except:
            e += 1
            print("Error: %d" % (i + 1))

    print("Checked %d reactions. %d reactions consist exception errors" % (i + 1, e))
    return 0


if __name__ == '__main__':
    main()
