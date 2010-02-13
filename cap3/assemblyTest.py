#!/usr/bin/env python
# encoding: utf-8
"""
assemblyTest.py

Created by Brant Faircloth on 2009-07-18.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import sqlite3
import assembly

def main():
    database = '/Volumes/Data/454_microsatellites/454_MSATS_7-3-09/test/Bog_copper/Blat/seq-seq-psl-blat.sqlite'
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    sequences = cur.execute('select sequence.id, sequence.seq_id, sequence.seq_seq from sequence').fetchall()
    # we have to read everything in because we cannot iterate as it locks
    # the database (allowing no further updates while iterating)
    count = 0
    for seq in sequences:
        a = assembly.Cluster(seq)
        # pass only the connection, we'll create a new cursor instance in the 
        # class to keep from overwriting the cursor instance (and its data)
        # here
        a.query(conn, cur)
        # if the base sequence returns data (i,e. is not masked/skipped)
        #pdb.set_trace()
        if a.data:
            a.filter()
            a.cap3()
            if a.cap3:
                a.baseContig('/Users/bcf/tmp/msatcommander-454/%s.contig.cap.ace' % a.base_name)
            if a.aligned:
                a.flag()
        count +=1
        print count
        #pdb.set_trace()


if __name__ == '__main__':
    main()