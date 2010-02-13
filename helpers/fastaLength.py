#!/usr/bin/env python
# encoding: utf-8
"""
fastaLength.py

Created by Brant Faircloth on 2009-09-23.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import re, time, pdb

def fLength(input):
    handle = open(input, 'rU')
    pdb.set_trace()
    lines = handle.read().count('>')
    #fasta = re.match('^$', handle.read())
    handle.close()
    print lines


if __name__ == '__main__':
    start = time.time()
    fLength('../sequence/Schisto_run1_FX5ZTWB03.sff.fna')
    end = time.time()
    print end - start