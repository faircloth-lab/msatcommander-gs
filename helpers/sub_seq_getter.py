#!/usr/bin/env python
# encoding: utf-8
"""
sub_seq_getter.py

Created by Brant Faircloth on 2009-12-17.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os, sys, pdb, time, MySQLdb, cPickle, optparse, ConfigParser
from Bio.SeqIO import QualityIO 

def interface():
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')
    p.add_option('--species', '-s', dest = 'species', action='store', \
type='string', default = None, help='The species for which to grab sequence')
    p.add_option('--trimmed', '-t', dest = 'trimmed', action='store', \
type='string', default = 'False', help='Get trimmed sequence')
    p.add_option('--msat', '-m', dest = 'msat', action='store', \
type='string', default = 'False', help='Get trimmed sequence')
    p.add_option('--all', '-a', dest = 'all', action='store', \
type='string', default = 'False', help='Get all sequences')
    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

def not_trimmed(cur, conf, options, sequence, qual):
    cur.execute('SELECT name FROM sequence WHERE cluster = %s', (options.species))
    data = cur.fetchall()
    dataset = set()
    for d in data:
        dataset.add(d[0])
    seqs = QualityIO.PairedFastaQualIterator(
    open(conf.get('Input','sequence'), "rU"), 
    open(conf.get('Input','qual'), "rU"))
    try:
        while seqs:
            record = seqs.next()
            if record.name in dataset:
                sequence.write('%s' % record.format('fasta'))
                qual.write('%s' % record.format('qual'))
    except StopIteration:
        pass
    qual.close()
    sequence.close()

def trimmed(cur, conf, species, msat, sequence, qual):
    #pdb.set_trace()
    if msat == 'True':
        cur.execute('''SELECT record FROM sequence where cluster = %s and msat = 1''', (species))
    elif msat == 'False':
        cur.execute('''SELECT record FROM sequence where cluster = %s''', (species))
    for p_rec in cur.fetchall():
        record = cPickle.loads(p_rec[0])
        if len(record) > 0:
            sequence.write('%s' % record.format('fasta'))
            qual.write('%s' % record.format('qual'))
    qual.close()
    sequence.close()

def main():
    start_time = time.time()
    options, arg = interface()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    if options.all == 'True' and options.trimmed == 'True':
        # get all the cluster names from the dbase
        cur.execute('SELECT distinct(cluster) from sequence')
        clusters = cur.fetchall()
        for c in clusters:
            if c[0] != None:
                #pdb.set_trace()
                out = c[0] + '_PARSED.fsa'
                q = c[0] + '_PARSED.qual'
                sequence = open(out,'w')
                qual = open(q,'w')
                trimmed(cur, conf, c[0], options.msat, sequence, qual)
    elif options.trimmed == 'False':
        out = options.species + '_PARSED.fsa'
        q = options.species + '_PARSED.qual'
        sequence = open(out,'w')
        qual = open(q,'w')
        not_trimmed(cur, conf, options.species, options.msat, sequence, qual)
    elif options.trimmed == 'True':
        out = options.species + '_PARSED.fsa'
        q = options.species + '_PARSED.qual'
        sequence = open(out,'w')
        qual = open(q,'w')
        trimmed(cur, conf, options, sequence, qual)
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'


if __name__ == '__main__':
    main()

