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
from math import log
from Bio import SeqIO

def db_write(c, db_row):
    c.execute("INSERT INTO blat (id, q_name, t_name, strand, \
    percent, score, match, mismatch, rep_match, \
    ns, q_num_insert, q_gap_bases, t_num_insert, t_gap_bases, q_size, q_start,\
    q_end, t_size, t_start, t_end) VALUES (:pkey, :qName, :tName, :strand, \
    :percent, :score, :match, :misMatch, :repMatch, :Ns, :qNumInsert, \
    :qGapBases, :tNumInsert, :tGapBases, :qSize, :qStart, :qEnd, :tSize, \
    :tStart, :tEnd)", db_row)
    # go ahead and commit since there will be very, very many writes
    c.commit()

def primaryRecord(c, pkey, seq_name, seq_seq, seq_len):
    c.execute("INSERT INTO sequence (id, seq_id, seq_seq, seq_len) VALUES (?,?,?,?)", (pkey, seq_name, seq_seq, seq_len))
    c.commit()

def primaryMatchCount(c, pkey, match_count):
    c.execute("UPDATE sequence SET match_count = ? WHERE id = ?", (match_count, pkey))
    c.commit()

def pslPercentId(psl, protein=False, mRNA=True):
    '''convert the percent sequence ID of an alignment from a single line of a 
    parsed PSL file.  Code adapted from 
    http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    millibad = 0
    if protein:
        sizeMul = 3
    else:
        sizeMul = 1
    qAliSize = sizeMul * (psl['qEnd'] - psl['qStart'])
    tAliSize = psl['tEnd'] - psl['tStart']
    aliSize = min(qAliSize, tAliSize)
    if aliSize <= 0:return 0
    else:
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
            if mRNA:
                sizeDif = 0;
            else:
                sizeDif = -sizeDif
        insertFactor = psl['qNumInsert']
        if not mRNA:
            insertFactor += psl['tNumInsert']
        total = (sizeMul * (psl['match'] + psl['repMatch'] + psl['misMatch']))
        if total != 0:
            milliBad = (1000 * (psl['misMatch']*sizeMul + insertFactor + round(3*log(1+sizeDif)))) / total
        percent = round(100 - milliBad * 0.1,0)
        psl['percent'] = percent
        return psl

def pslScore(psl, sizeMul = 1):
    '''convert the score of an alignment from a single line of a 
    parsed PSL file.  Code adapted from 
    http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    score = sizeMul * (psl['match'] + (psl['repMatch'] >> 1)) - sizeMul * psl['misMatch'] - psl['qNumInsert'] - psl['tNumInsert']
    psl['score'] = score
    return psl

def no_hits(c, pkey, seq_name):
    #fieldnames = ["subject_name", "match_percent", "match_length", "mismatches", "gap_openings", "query_start", "query_end", "subject_start", "subject_end", "e_score", "bit_score"]
    fieldnames = ['match', 'misMatch', 'repMatch', 'Ns', 'qNumInsert', 'qGapBases', 'tNumInsert', 'tGapBases', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'qStarts', 'tStarts']
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
    # score of 45 roughyly approximates e_score ~ 1e-15
    # score of 50 roughly approximates e_score ~ 1e-20
    score = 45
    fieldnames = ['match', 'misMatch', 'repMatch', 'Ns', 'qNumInsert', 'qGapBases', 'tNumInsert', 'tGapBases', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes','qStarts', 'tStarts']
    match_count = 0
    self_match = False
    for row in blat_result:
        row = row.split()
        # merge fieldnames with row data
        db_row = dict(zip(fieldnames, row))
        # covert strings to ints
        for f in ['match', 'misMatch', 'repMatch', 'Ns', 'qNumInsert', 'qGapBases', 'tNumInsert', 'tGapBases', 'qSize', 'qStart', 'qEnd', 'tSize', 'tStart', 'tEnd']:
            db_row[f] = int(db_row[f])
        # compute score
        db_row = pslScore(db_row)
        # compute percent id
        db_row = pslPercentId(db_row)
        #pdb.set_trace()
        # if this is NOT a self-match
        if not (db_row['qName'] == db_row['tName']) and db_row['score'] >= score:
            # update the dependent table with the foreign key
            db_row['pkey'] = pkey
            db_write(c, db_row)
            match_count += 1
        elif (db_row['qName'] == db_row['tName']):
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
        #blat_result, blat_error = subprocess.Popen('/Users/bcf/bin/i386/gfclient -t=DNA -q=DNA -minScore=0 -minIdentity=0 -out=psl -nohead 127.0.0.1 8888 / stdin stdout', shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE).communicate(input=fasta_record)
        #pdb.set_trace()
        blat_result, blat_error = subprocess.Popen('/Users/bcf/bin/i386/blat /Users/bcf/Data/454_msat/blat/Bog_Copper.softmask.clean.2bit stdin -mask=lower -out=psl -noHead stdout', shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE).communicate(input=fasta_record)
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
    q_name text,
    t_name text,
    strand text,
    percent real,
    score real,
    match int,
    mismatch int,
    rep_match int,
    ns int,
    q_num_insert int,
    q_gap_bases int,
    t_num_insert int,
    t_gap_bases int,
    q_size int,
    q_start int,
    q_end int,
    t_size int,
    t_start int,
    t_end int,
    foreign key (id) references sequence(id)
    )''')

def main():
    # start server - don't start alignment in repeat region (not sure if this is working correctly)
    # /Users/bcf/bin/i386/gfserver start 127.0.0.1 8888 -mask /Users/bcf/data/genome/mm9/mm9.2bit
    # start server - more verbose and more spurious matches with -stepSize=5
    # /Users/bcf/bin/i386/gfserver start 127.0.0.1 8888 -stepSize=5 /Users/bcf/data/genome/mm9/mm9.2bit
    # manual gfServer command
    # /Users/bcf/bin/i386/gfclient -t=DNA -q=DNA -minScore=0 -minIdentity=0 -out=psl 127.0.0.1 8888 / ~/tmp/tmp.fa stdout
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    input_fasta = '/Users/bcf/data/454_msat/blat/Bog_Copper_2EE_454_MSATs.fa.clean'
    c = sqlite3.connect('/Users/bcf/data/454_msat/blat/seq-seq-psl-blat.sqlite')
    addTable(c)
    blast_screen(c, input_fasta)
    c.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'


if __name__ == '__main__':
    main()

