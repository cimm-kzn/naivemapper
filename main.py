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
import argparse
import sys
import traceback
from CGRtools.RDFread import RDFread
from CGRtools.RDFwrite import RDFwrite
from core.fragger import Fragger


def mapper_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = RDFwrite(kwargs['output'])

    fragger = Fragger(_min=1, _max=6)
    err = 0
    num = 0

    for num, data in enumerate(inputdata.readdata(), start=1):
        if kwargs['debug'] or num % 10 == 1:
            print("reaction: %d" % num, file=sys.stderr)
        try:
            a = some_magic
            outputdata.writedata(a)
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)

    print('%d from %d reactions processed' % (num - err, num), file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(description="", epilog="", prog='naivemapper')
    parser.add_argument("--version", "-v", action="version", version="0.1", default=False)
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'),
                        help="RDF outputfile")
    parser.add_argument("--model", "-m", type=str, default='model.dat',
                        help="path to model save|restore file")
    parser.add_argument("--fit", "-f", action='store_true', help="Fit model")

    parser.add_argument("--debug", action='store_true', help="debug mod")
    return parser


def main():
    parser = parse_args()
    args = parser.parse_args()
    mapper_core(**vars(args))

if __name__ == '__main__':
    main()
