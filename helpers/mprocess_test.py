#!/usr/bin/env python
# encoding: utf-8
"""
mprocess_test.py

Created by Brant Faircloth on 2009-08-24.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os
import pdb
import time
import ConfigParser
import multiprocessing
from Bio.SeqIO import QualityIO
from Bio import pairwise2


def worker(record, count):
    """worker function"""
    seq, tag = str(record.seq), str(record.seq)
    seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
    #name = multiprocessing.current_process().name
    #print 'Worker', name, str(record.seq)
    print "Parent: ", os.getppid(), "Child: ", os.getpid(), "Count: ", count
    return

if __name__ == '__main__':
    start_time = time.time()
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    #jobs = []
    record = QualityIO.PairedFastaQualIterator(open(conf.get('Input','sequence'), "rU"), open(conf.get('Input','qual'), "rU"))
    mproc=True
    if mproc == True:
        count = 0
        try:
            while count < 500:
                #pdb.set_trace()
                jobs = []
		for i in range(multiprocessing.cpu_count()):
                    count +=1
                    p = multiprocessing.Process(target=worker, args=(record.next(),count))
                    jobs.append(p)
                    p.start()
                #p.join()
        except StopIteration:
            pass
    else:
        count = 0
	try:
            while count < 500:
                count +=1
                worker(record.next(), count)
        except StopIteration:
            pass
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

