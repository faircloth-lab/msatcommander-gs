#!/usr/bin/env python
# encoding: utf-8
"""
subclassAce.py

Created by Brant Faircloth on 2009-07-16.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
from Bio.Sequencing import Ace


def write(input, output, wrap = 60):
    '''write(input, output)
    
    where input is a Contig object and output is a location on the directory 
    tree.
    
    wrap is an optional parameter used to determine the line-wrapping width of
    sequence and quality data.  Default = 60.
    
    This function writes a Contig object to an ACE formmatted file at `output`
    '''
    h = open(output, 'w')
    h.write('AS %s %s\n\n' % (input.ncontigs, input.nreads))
    for c in input.contigs:
        seq_wrapped = [c.sequence[i:i+wrap] for i in xrange(0, c.nbases, wrap)]
        h.write('CO %s %s %s %s %s\n%s\n\nBQ\n' % (c.name, c.nbases, c.nreads,\
        c.nsegments, c.uorc, '\n'.join(seq_wrapped)))
        for i in xrange(0,len(c.quality),wrap/3):
            h.write('%s\n' % ' '.join(['%s' % q for q in c.quality[i:i+wrap/3]]))
        h.write('\n')
        for a in c.af:
            h.write('AF %s %s %s\n' % (a.name, a.coru, a.padded_start))
        for b in c.bs:
            h.write('BS %s %s %s\n' % (b.padded_start, b.padded_end, b.name))
        h.write('\n')
        for r in c.reads:
            seq_wrapped = [r.rd.sequence[i:i+wrap] for i in xrange(0, \
            len(r.rd.sequence), wrap)]
            h.write('RD %s %s %s %s\n%s\n\n' % (r.rd.name, r.rd.padded_bases, \
            r.rd.info_items, r.rd.read_tags, '\n'.join(seq_wrapped)))
            h.write('QA %s %s %s %s\n' % (r.qa.qual_clipping_start, \
            r.qa.qual_clipping_end, r.qa.align_clipping_start, r.qa.align_clipping_end))
            if r.ds.chromat_file == '':
                h.write('DS \n\n')
            else:
                h.write('DS CHROMAT_FILE: %s PHD_FILE: %s TIME: %s CHEM: %s \
                DYE: %s TEMPLATE: %s DIRECTION: %s\n\n' % (r.ds.chromat_file, \
                r.ds.phd_file, r.ds.time, r.ds.chem, r.ds.dye, r.ds.template, \
                r.ds.direction))
    h.close()


def main():
    base_name = 'FX5ZTWB02D1DFX'
    #contigs = Ace.parse(open('/Users/bcf/Tmp/tmp2.fa.cap.ace'))
    c = Ace.read(open('/Users/bcf/Tmp/tmp2.fa.cap.ace'))
    '''for c in contigs:
        for r in c.reads:
            if r.rd.name == base_name:
                contig = c
                break
            else:
                pass'''
    write(c, '/Users/bcf/Tmp/tmp2_rewrite.fa.cap.ace')
    pdb.set_trace()

if __name__ == '__main__':
    main()

