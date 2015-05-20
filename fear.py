#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of FEAR (Find Errors in Automapped Reactions).
#
#  FEAR is free software; you can redistribute it and/or modify
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
import Tkinter as tk
import os

from tkFileDialog import askopenfilename, asksaveasfilename

from FEAR.CGR import CGR
from FEAR.RDFread import RDFread
from FEAR.version import version
from FEAR.RDFwrite import RDFwrite


class GUI:
    def __init__(self, master):
        self.files = dict()

        self.last_used_directory = '.'
        frame = tk.Frame(master)
        frame.pack()

        self.button = tk.Button(frame, text="input RDF", command=self.ifilename)
        self.button.pack(side=tk.LEFT)

        self.hi_there = tk.Button(frame, text="output RDF", command=self.ofilename)
        self.hi_there.pack(side=tk.LEFT)

        self.button = tk.Button(frame, text="QUIT", command=frame.quit)
        self.button.pack(side=tk.LEFT)

    def ifilename(self):
        self.files['input'] = askopenfilename(initialdir=self.last_used_directory)

    def ofilename(self):
        fullfilename = asksaveasfilename(initialdir=self.last_used_directory)
        if os.path.splitext(fullfilename)[1].lower() == ".rdf":
            self.files['output'] = fullfilename
        else:
            self.files['output'] = '.'.join([fullfilename, 'rdf'])


def main():
    rawopts = argparse.ArgumentParser(description="Find Errors in Automapped Reactions",
                                      epilog="Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default='input.rdf', help="input RDF ")
    rawopts.add_argument("--gui", action='store_true', help="start GUI")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="output RDF")
    options = vars(rawopts.parse_args())

    config = os.path.join(os.path.expanduser('~'), '.fear', 'conf')
    if os.path.exists(config):
        print ('ggg')
    else:
        print ('fff')

    if options['gui']:
        root = tk.Tk()
        app = GUI(root)
        root.mainloop()
        options.update(app.files)
        if options['output'] == 'output.rdf':
            path, filename = os.path.split(options['input'])
            options['output'] = os.path.join(path, 'output.rdf')
        print (options)
        return 0

    rdf = RDFread(options['input'])
    if not rdf.chkRDF():
        print('rdf incorrect')
        return 0

    corr = RDFwrite(options['output'])
    incorr = RDFwrite(options['output'][:-3]+'err.rdf')
    result = []
    stats = []
    calc = CGR()
    e = 0
    for i, data in enumerate(rdf.readdata()):
        if i % 1000 == 0 and i:
            print ("reaction: %d" % (i + 1))
        #res = calc.firstcgr(data)
        try:
            res = calc.firstcgr(data, i +1)
            if res:
                corr.write(data)
            else:
                incorr.write(data)
        except:
            e += 1
            print ("Error: %d" % (i+1))

    print("Checked %d reactions. %d reactions consist exception errors" % (i+1, e))
    return 0


if __name__ == '__main__':
    main()
