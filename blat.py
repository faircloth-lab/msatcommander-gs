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
import sys, pdb, time, MySQLdb, subprocess, ConfigParser, tempfile, cPickle, multiprocessing, textwrap
from math import log
#from Bio import SeqIO

def db_write(cur, db_row):
    cur.execute('''INSERT INTO blat (id, q_name, t_name, strand, percent, 
    score, matches, mismatches, rep_match, ns, q_num_insert, q_gap_bases, 
    t_num_insert, t_gap_bases, q_size, q_start, q_end, t_size, t_start, 
    t_end) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
    (db_row['pkey'], db_row['qName'], db_row['tName'], db_row['strand'], 
    db_row['percent'], db_row['score'], db_row['match'], db_row['misMatch'], 
    db_row['repMatch'], db_row['Ns'], db_row['qNumInsert'], 
    db_row['qGapBases'], db_row['tNumInsert'], db_row['tGapBases'], 
    db_row['qSize'], db_row['qStart'], db_row['qEnd'], db_row['tSize'], 
    db_row['tStart'], db_row['tEnd']))
    # go ahead and commit since there will be very, very many writes
    #cur.commit()

def primaryMatchCount(cur, pkey, match_count):
    cur.execute('''UPDATE sequence_test SET match_count = %s WHERE id = %s''', (match_count, pkey))

def pslPercentId(psl, protein=False, mRNA=True):
    '''
    convert the percent sequence ID of an alignment from a single line of a 
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
    '''
    convert the score of an alignment from a single line of a parsed PSL file.
    Code adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    score = sizeMul * (psl['match'] + (psl['repMatch'] >> 1)) - sizeMul * psl['misMatch'] - psl['qNumInsert'] - psl['tNumInsert']
    psl['score'] = score
    return psl

def parse_blat(cur, pkey, blat_result, seq_name, seq_seq):
    # score of 40 roughly approximates e_score ~ 1e-10
    # score of 45 roughly approximates e_score ~ 1e-15
    # score of 50 roughly approximates e_score ~ 1e-20
    #pdb.set_trace()
    score = 40
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
        # if this is NOT a self-match
        if not (db_row['qName'] == db_row['tName']) and db_row['score'] >= score:
            # update the dependent table
            db_row['pkey'] = pkey
            db_write(cur, db_row)
            match_count += 1
        elif (db_row['qName'] == db_row['tName']):
            # self match
            pass
    # update the primary table with the match count
    primaryMatchCount(cur,pkey,match_count)
    #pdb.set_trace()

def worker(tb, pkey, name, seq_trimmed):
        conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
        db="454_msatcommander")
        cur = conn.cursor()
        sequence = ('>%s\n%s\n' % (name, seq_trimmed))
        #pdb.set_trace()
        blat_result, blat_error = subprocess.Popen('/Users/bcf/bin/i386/blat %s stdin -mask=lower -out=psl -noHead stdout' % tb, shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate(input=sequence)
        if blat_error:
            pdb.set_trace()
        if blat_result and not blat_error:
            # drop the ending newline
            blat_result = blat_result.split('\n')[:-1]
            parse_blat(cur, pkey, blat_result, name, seq_trimmed)
            #print 'parsing'
        elif not blat_result and not blat_error:
            #pdb.set_trace()
            primaryMatchCount(cur, pkey, 0)
        elif blat_error:
            print "blat processing error: %s" % blat_error
            #pdb.set_trace()
        conn.commit()
        cur.close()
        conn.close()

def createBlatTable(cur):
    '''add a table to the database'''
    try:
        # if previous tables exist, drop them
        # TODO: fix createDbase() to drop tables safely
        cur.execute('''DROP TABLE blat''')
    except:
        pass
    cur.execute('''CREATE TABLE blat (id INT UNSIGNED NOT NULL, q_name 
        VARCHAR(100), t_name VARCHAR(100), strand VARCHAR(1), percent 
        DECIMAL(4,1), score DECIMAL(4,1), matches SMALLINT UNSIGNED, mismatches 
        SMALLINT UNSIGNED, rep_match SMALLINT UNSIGNED, ns SMALLINT UNSIGNED, 
        q_num_insert SMALLINT UNSIGNED, q_gap_bases SMALLINT UNSIGNED, 
        t_num_insert SMALLINT UNSIGNED, t_gap_bases SMALLINT UNSIGNED, q_size 
        SMALLINT UNSIGNED, q_start SMALLINT UNSIGNED, q_end SMALLINT 
        UNSIGNED, t_size SMALLINT UNSIGNED, t_start SMALLINT UNSIGNED, t_end 
        SMALLINT UNSIGNED, INDEX blat_id (id))''')
        
def updateSequenceTable(cur):
    #pdb.set_trace()
    try:
        cur.execute('''ALTER TABLE sequence_test ADD COLUMN match_count SMALLINT UNSIGNED''')
        print 'Added match_count column'
    except MySQLdb._mysql.OperationalError, e:
        if e[0] == 1060:
            cur.execute('''UPDATE sequence_test SET match_count = NULL''')
            print 'Zeroed-out match_count column'
    except:
        print "Error adding match_count to sequence table"
        sys.exit()


def fasta(name, seq):
    return '>%s\n%s\n' % (name, seq)

def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    # create a table for the output
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
    db="454_msatcommander")
    cur = conn.cursor()
    createBlatTable(cur)
    updateSequenceTable(cur)
    # get clusters from conf file
    cluster = [c[1] for c in conf.items('Clusters')]
    # create a temp directory for results
    tdir = tempfile.mkdtemp(prefix='msatcommander-454-')
    #print tdir
    # get sequences from Dbase by cluster
    for c in cluster:
        # get all sequence to build the 2bit file
        cur.execute('''SELECT id, name, seq_trimmed FROM sequence_test WHERE cluster = %s''', (c,))
        sequences = cur.fetchall()
        if sequences:
            tf = tempfile.mkstemp(prefix='clean-%s-' % c, suffix='.fa', dir=tdir)
            tf_handle = open(tf[1], 'w')
            for s in sequences:
                tf_handle.write(fasta(s[1],s[2]))
            tf_handle.close()
            # for each cluster generate the 2bit file for blat
            tb = tf[1] + '.2bit'
            #tb_handle = open(tb,'w')
            print "Building %s twobit file" % c
            #two_bit = subprocess.Popen('/Users/bcf/bin/i386/faToTwoBit %s stdout' % tf[1], shell=True, stdout=tb_handle, stderr = subprocess.PIPE).communicate()
            two_bit = subprocess.Popen('/Users/bcf/bin/i386/faToTwoBit %s %s' % (tf[1], tb), shell=True, stderr = subprocess.PIPE).communicate()
            #tb_handle.close()
            # get only msat-containing sequences
            cur.execute('''SELECT id, name, seq_trimmed FROM sequence_test WHERE cluster = %s and msat = 1''', (c,))
            sequences = cur.fetchall()
            if conf.getboolean('Multiprocessing', 'MULTIPROCESSING'):
                # get num processors
                n_procs = conf.get('Multiprocessing','processors')
                if n_procs == 'Auto':
                    n_procs = multiprocessing.cpu_count() - 1
                else:
                    n_procs = int(n_procs)
                print 'Multiprocessing.  Number of processors = ', n_procs
                threads = []
                index = 0
                while index < len(sequences):
                #while index < 50:
                    if len(threads) < n_procs:
                        p = multiprocessing.Process(target=worker, args=(tb, sequences[index][0], sequences[index][1], sequences[index][2]))
                        p.start()
                        threads.append(p)
                        index += 1
                        #if index%1000 == 0:
                        #    print "Sequence # %s, %s" % (index, time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(time.time())))
                    else:
                        for t in threads:
                            if not t.is_alive():
                                threads.remove(t)
                # clean up any running processes before moving along to next cluster
                while threads:
                    for t in threads:
                        if not t.is_alive():
                            threads.remove(t)
            else:
                print 'Not using multiprocessing'
                index = 0
                #while index < len(sequences):
                while index < 50:
                    # send 0th item of list to worker and remove
                    worker(tb, sequences[index][0], sequences[index][1], sequences[index][2])
                    index += 1
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'
    #pdb.set_trace()


if __name__ == '__main__':
    main()

