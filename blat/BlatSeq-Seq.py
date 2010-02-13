#!/usr/bin/env python
# encoding: utf-8
"""
BlatGenomicScreener.py

Created by Brant Faircloth on 2009-05-22.
Copyright (c) 2009 Brant Faircloth. All rights reserved.

Takes fasta file as input, and screens content against genomic data
(.2bit file) using Blat as the sequence comparison engine and writing data to 
a SQLite dbase.

"""
import sys
import pdb
import time
import sqlite3
import subprocess
from Bio import SeqIO

def db_write(c, db_row):
    c.execute("INSERT INTO blat (id, query_name, subject_name, match_percent, \
    match_length, mismatches, gap_openings, query_start, query_end, \
    subject_start, subject_end, e_score, bit_score) VALUES \
    (:pkey, :query_name, :subject_name, :match_percent, :match_length, :mismatches, \
    :gap_openings, :query_start, :query_end, :subject_start, :subject_end, \
    :e_score, :bit_score)", db_row)
    # go ahead and commit since there will be very, very many writes
    c.commit()

def primaryRecord(c, pkey, seq_name, seq_seq, seq_len):
    c.execute("INSERT INTO sequence (id, seq_id, seq_seq, seq_len) VALUES (?,?,?,?)", (pkey, seq_name, seq_seq, seq_len))
    c.commit()

def primaryMatchCount(c, pkey, match_count):
    c.execute("UPDATE sequence SET match_count = ? WHERE id = ?", (match_count, pkey))
    c.commit()
    
def no_hits(c, pkey, seq_name):
    fieldnames = ["subject_name", "match_percent", "match_length", "mismatches", "gap_openings", "query_start", "query_end", "subject_start", "subject_end", "e_score", "bit_score"]
    seq_name_split = seq_name.split('_')
    query_name_ucsc = seq_name
    row = [None]*len(fieldnames)
    db_row = dict(zip(fieldnames,row))
    db_row['query_name'],db_row['query_ucsc'],db_row['subject_ucsc'] = seq_name, query_name_ucsc, None
    db_write(c, db_row)

def parse_blat(c, pkey, blat_result, seq_name, seq_seq):
    # make a primary entry for the sequence
    primaryRecord(c, pkey, seq_name, str(seq_seq), len(seq_seq))
    # proceed with the parsing
    e_score = 1e-15
    fieldnames = ["query_name", 
    "subject_name", 
    "match_percent", 
    "match_length", 
    "mismatches", 
    "gap_openings", 
    "query_start", 
    "query_end", 
    "subject_start", 
    "subject_end", 
    "e_score", 
    "bit_score"
    ]
    match_count = 0
    self_match = False
    for row in blat_result:
        row = row.split()
        # merge fieldnames with row data
        db_row = dict(zip(fieldnames, row))
        # if this is NOT a self-match
        if not (db_row["query_name"] == db_row["subject_name"]) and float(db_row['e_score']) <= e_score:
            # update the dependent table with the foreign key
            db_row["pkey"] = pkey
            db_write(c, db_row)
            match_count += 1
        elif (db_row["query_name"] == db_row["subject_name"]):
            # self match
            pass
    # update the primary table with the match count
    primaryMatchCount(c,pkey,match_count)
    #pdb.set_trace()

def blast_screen(c, input_fasta):
    pkey = 1
    for seq in SeqIO.parse(open(input_fasta,'rU'),'fasta'):
        # this is slightly redundant, but it keeps us from having to deal with
        # line endings and whatnot (if there are line ending issues).
        fasta_record = seq.format("fasta")
        # os.popen() is deprecated
        blat_result, blat_error = subprocess.Popen('/Users/bcf/bin/i386/gfclient -t=DNA -q=DNA -minScore=0 -minIdentity=0 -out=blast8 127.0.0.1 8888 / stdin stdout', shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE).communicate(input=fasta_record)
        #pdb.set_trace()
        if blat_result and not blat_error:
            # drop the ending newline
            blat_result = blat_result.split('\n')[:-1]
            parse_blat(c, pkey, blat_result, seq.name, seq.seq)
        elif not blat_result and not blat_error:
            # this shouldn't happen any longer given no 'N' in seq.seq
            #pdb.set_trace()
            primaryMatchCount(c, pkey, None)
        elif blat_error:
            print "blat processing error: %s" % blat_error
        pkey += 1
        if pkey%1000 == 0:
            print (("Sequence # %s, %s") % (pkey, time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(time.time()))))

def addTable(c):
    '''add a table to the database'''
    try:
        # if previous tables exist, drop them
        # TODO: fix createDbase() to drop tables safely
        c.execute('''drop table blat''')
        # commit the change
        c.commit()
    except:
        pass
    # create the primers results table
    c.execute('''create table sequence (
    id integer primary key,
    seq_id text,
    seq_len int,
    seq_seq text,
    seq_qual text, 
    seq_file_path text,
    match_count int
    )''')
    c.execute('''create table blat (
    id integer,
    query_name text,
    subject_name text,
    match_percent real,
    match_length int,
    mismatches int,
    gap_openings int,
    query_start int,
    query_end int,
    subject_start int,
    subject_end int,
    e_score real,
    bit_score real,
    foreign key (id) references sequence(id)
    )''')

def main():
    # start server - don't start alignment in repeat region (not sure if this is working correctly)
    # /Users/bcf/bin/i386/gfserver start 127.0.0.1 8888 -mask /Users/bcf/data/genome/mm9/mm9.2bit
    # start server - more verbose and more spurious matches with -stepSize=5
    # /Users/bcf/bin/i386/gfserver start 127.0.0.1 8888 -stepSize=5 /Users/bcf/data/genome/mm9/mm9.2bit
    # manual gfServer command
    # /Users/bcf/bin/i386/gfclient -t=DNA -q=DNA -minScore=0 -minIdentity=0 -out=blast8 127.0.0.1 8888 / tmp.fa stdout
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    input_fasta = '/Users/bcf/data/454_msat/MUS_454_MSATs.fa.clean'
    c = sqlite3.connect('seq-genome.sqlite')
    addTable(c)
    blast_screen(c, input_fasta)
    c.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'


if __name__ == '__main__':
    main()

