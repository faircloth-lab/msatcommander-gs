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
import os, sys, pdb, time, MySQLdb, subprocess, ConfigParser, tempfile, cPickle, multiprocessing, textwrap, optparse
from math import log
#from Bio import SeqIO

def db_write(cur, db_row):
    cur.execute('''INSERT INTO blat (id, t_id, strand, percent, 
    length, mismatches, q_gaps, q_size, q_start, q_end, t_size, t_start, 
    t_end, e_score, bit_score) 
    VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
    (db_row['qName'], db_row['tName'], db_row['strand'], 
    db_row['percent'], db_row['length'], db_row['misMatch'], db_row['qGaps'], 
    db_row['qSize'], db_row['qStart'], db_row['qEnd'], db_row['tSize'], 
    db_row['tStart'], db_row['tEnd'], db_row['eScore'], db_row['bitScore']))
    # go ahead and commit since there will be very, very many writes
    #cur.commit()

def primaryMatchCount(cur, pkey, match_count):
    cur.execute('''UPDATE sequence SET match_count = %s WHERE id = %s''', (match_count, pkey))

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


def parse_blat(conn, cur, pkey, blat_result, seq_name, seq_seq):
    percent = 90.0
    length = 25
    e_score = 1e-10
    fieldnames = ['qName','tName', 'percent', 'length', 'misMatch', 'qGaps', 'qStart', 'qEnd', 'tStart', 'tEnd', 'eScore', 'bitScore']
    match_count = 0
    self_match = False
    unique_matches = []
    for row in blat_result:
        row = row.split()
        # merge fieldnames with row data
        db_row = dict(zip(fieldnames, row))
        # covert strings to ints
        for f in ['length', 'misMatch', 'qGaps', 'qStart', 'qEnd', 'tStart', 'tEnd',]:
            db_row[f] = int(db_row[f])
        for f in ['percent', 'eScore', 'bitScore']:
            db_row[f] = float(db_row[f])
        # determine strand
        if db_row['tStart'] > db_row['tEnd']:
            db_row['tStart'], db_row['tEnd'] = db_row['tEnd'], db_row['tStart']
            db_row['strand'] = '-'
        else:
            db_row['strand'] = '+'
        # fix/add some other metrics
        db_row['percent'] = round(db_row['percent'],1)
        db_row['qSize'] = db_row['qEnd'] - db_row['qStart']
        db_row['tSize'] = db_row['tEnd'] - db_row['tStart']
        # if this is NOT a self-match
        if not (db_row['qName'] == db_row['tName']) and db_row['percent'] >= percent and db_row['length'] >= length and db_row['eScore'] <= e_score:
            # update the dependent table
            db_row['pkey'] = pkey
            db_write(cur, db_row)
            # it appears we need to commit here to avoid a mysterious deadlock issue
            conn.commit()
            if db_row['tName'] not in unique_matches:
                unique_matches.append(db_row['tName'])
                match_count += 1
        elif (db_row['qName'] == db_row['tName']):
            # self match
            pass
    # update the primary table with the match count
    primaryMatchCount(cur,pkey,match_count)
    # it appears we need to commit here to avoid a mysterious deadlock issue
    conn.commit()
    #pdb.set_trace()

def worker(conf, tb, pkey, name, seq_trimmed):
        conn = MySQLdb.connect(user=conf.get('Database','USER'), 
            passwd=conf.get('Database','PASSWORD'), 
            db=conf.get('Database','DATABASE'))   
        cur = conn.cursor()
        # DONE: change to pkey right here
        sequence = ('>%s\n%s\n' % (pkey, seq_trimmed))
        #pdb.set_trace()
        blat_result, blat_error = subprocess.Popen('/Users/bcf/bin/i386/blat %s stdin -mask=lower -out=blast8 -noHead stdout' % tb, shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate(input=sequence)
        if blat_error:
            pdb.set_trace()
        if blat_result and not blat_error:
            # drop the ending newline
            blat_result = blat_result.split('\n')[:-1]
            parse_blat(conn, cur, pkey, blat_result, name, seq_trimmed)
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
        # DONE:  change q_name and t_name columns to q_id and t_id when ref changes
        # DONE:  change referencing to foreign key type w/ InnoDB
    cur.execute('''CREATE TABLE blat (id INT UNSIGNED NOT NULL, t_id INT 
        UNSIGNED NOT NULL, strand VARCHAR(1), percent 
        DECIMAL(4,1), length SMALLINT unsigned, mismatches SMALLINT UNSIGNED, 
        q_gaps SMALLINT UNSIGNED, q_size SMALLINT UNSIGNED, q_start SMALLINT 
        UNSIGNED, q_end SMALLINT UNSIGNED, t_size SMALLINT UNSIGNED, t_start 
        SMALLINT UNSIGNED, t_end SMALLINT UNSIGNED, e_score FLOAT, bit_score 
        FLOAT, FOREIGN KEY (id) REFERENCES sequence (id), FOREIGN KEY 
        (t_id) REFERENCES sequence (id)) ENGINE=InnoDB''')
        
def updateSequenceTable(cur):
    #pdb.set_trace()
    try:
        cur.execute('''ALTER TABLE sequence ADD COLUMN match_count SMALLINT UNSIGNED''')
        cur.execute('''ALTER TABLE sequence ADD INDEX match_count_idx (match_count)''')
        print 'Added match_count column and index'
    except MySQLdb._mysql.OperationalError, e:
        if e[0] == 1060:
            cur.execute('''UPDATE sequence SET match_count = NULL''')
            print 'Zeroed-out match_count column'
    except:
        print "Error adding match_count to sequence table"
        sys.exit()


def fasta(name, seq):
    return '>%s\n%s\n' % (name, seq)

def motd():
    motd = '''
    ##############################################################
    #                     msatcommander 454                      #
    #                                                            #
    # - parsing and error correction for sequence tagged primers #
    # - microsatellite identification                            #
    # - sequence pooling                                         #
    # - primer design                                            #
    #                                                            #
    # Copyright (c) 2009 Brant C. Faircloth & Travis C. Glenn    #
    ##############################################################\n
    '''
    print motd

def interface():
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

def q_runner(n_procs, list_item, function, *args):
    '''generic function used to start worker processes'''
    task_queue      = multiprocessing.Queue()
    results_queue   = multiprocessing.JoinableQueue()
    if args:
        arguments = (task_queue, results_queue,) + args
    else:
        arguments = (task_queue, results_queue,)
    results = []
    # reduce processer count if proc count > files
    if len(list_item) < n_procs:
        n_procs = len(list_item)
    for l in list_item:
        task_queue.put(l)
    for _ in range(n_procs):
        p = multiprocessing.Process(target = function, args = arguments).start()
        print 'Starting %s' % function
    for _ in range(len(list_item)):
        # indicated done results processing
        results.append(results_queue.get()) 
        results_queue.task_done()
    #tell child processes to stop
    for _ in range(n_procs):
        task_queue.put('STOP')
    # join the queue until we're finished processing results
    results_queue.join()
    # not closing the Queues caused me untold heartache and suffering
    task_queue.close()
    results_queue.close()
    return results

def worker1(input, output, conf, tdir, msat=False):
    '''docstring for worker1'''
    for c in iter(input.get, 'STOP'):
        # get all sequences to build the 2bit file
        conn = MySQLdb.connect(user=conf.get('Database','USER'), 
            passwd=conf.get('Database','PASSWORD'), 
            db=conf.get('Database','DATABASE'))   
        cur = conn.cursor()
        if msat:
            pf = 'msat'
            cur.execute('''SELECT id, name, seq_trimmed FROM sequence WHERE cluster = %s and msat = 1''', (c,))
        else:
            pf = 'clean'
            cur.execute('''SELECT id, name, seq_trimmed FROM sequence WHERE cluster = %s''', (c,))
        sequences = cur.fetchall()
        if sequences:
            tf = tempfile.mkstemp(prefix='%s-%s-' % (pf,c), suffix='.fa', dir=tdir)
            tf_handle = open(tf[1], 'w')
            for s in sequences:
                # DONE:  change s[1] to s[0] right here to add target id
                tf_handle.write(fasta(s[0],s[2]))
            print 'closing file'
            tf_handle.close()
        cur.close()
        conn.close()
        output.put(tf)

def worker2(input, output):
    '''docstring for worker2'''
    for f in iter(input.get, 'STOP'):
        print '2 bit process started'
        # for each cluster generate the 2bit file for blat
        mf = f + '.masked'
        tb = f + '.2bit'
        print "Building twobit file from\t=>\t\t\t\t%s" % mf
        two_bit = subprocess.Popen('/Users/bcf/bin/i386/faToTwoBit %s %s' % (mf, tb), shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        if two_bit[1]:
            output.put(two_bit[1])
        else:
            output.put(tb)
        print '2 bit process started'


def worker3(input, output, tdir):
    '''docstring for worker3'''
    for f in iter(input.get, 'STOP'):
        print 'blat process started'
        db, query   = f
        pre     = os.path.basename(query).split('-')[1]
        out     = tempfile.mkstemp(prefix='blat-%s-' % pre, suffix='.blat', dir=tdir)
        blat_error = subprocess.Popen('/Users/bcf/bin/i386/blat %s %s -mask=lower -out=blast8 -noHead %s' % (db, query, out[1]), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate(None)
        if blat_error[1]:
            output.put(blat_error[1])
        else:
            output.put(out)

def worker4(input, output, conf, file):
    for chunk in iter(input.get, 'STOP'):
        # we're getting a file chunk to read
        handle = open(file)
        handle.seek(chunk[0])
        data = handle.readlines(chunk[1])
        handle.close()
        parseBlat(conf, data    )

def repeatmasker(n_procs, fastas):
    for f in fastas:
        if n_procs > 1:
            print "Masking low complexity regions with RepeatMasker (multicore) in\t=>\t%s" % f[1]
            lc_masker_out, lc_masker_err = subprocess.Popen('/Users/bcf/bin/RepeatMasker/repeatmasker -qq -noint -xsmall -pa %s %s' % (n_procs - 1, f[1]), shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        else:
            print "Masking low complexity regions in\t=>\t%s with RepeatMasker" % f[1]
            lc_masker_out, lc_masker_err = subprocess.Popen('/Users/bcf/bin/RepeatMasker/repeatmasker -qq -noint -xsmall %s' % f[1], shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)

def chunker(file, size=1024*1024*5):
    '''chunking code borrowed and slightly modified from
    http://effbot.org/zone/wide-finder.htm#a-multi-threaded-python-solution'''
    handle = open(file)
    while 1:
        start = handle.tell()
        handle.seek(size,1)
        s = handle.readline()
        yield start, handle.tell() - start
        if not s:
            break
    handle.close()

def main():
    start_time = time.time()
    options, arg = interface()
    motd()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    # TODO:  increase number of processes because each is relatively
    # lightweight?
    if conf.getboolean('Multiprocessing', 'MULTIPROCESSING'):
        # get num processors
        n_procs = conf.get('Multiprocessing','processors')
        if n_procs == 'Auto':
            n_procs = multiprocessing.cpu_count() - 1
        else:
            n_procs = int(n_procs)
        print 'Multiprocessing.  Number of processors = ', n_procs
    else:
        n_procs = 1
        print 'Not using multiprocessing. Number of processors = ', n_procs
    # create a table for the output
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    createBlatTable(cur)
    updateSequenceTable(cur)
    conn.commit()
    # get clusters from the dbase
    cur.execute('''SELECT distinct(cluster) FROM SEQUENCE''')
    c_data = cur.fetchall()
    cluster = [c[0] for c in c_data]
    # create a temp directory for results
    tdir = tempfile.mkdtemp(prefix='msatcommander-454-')
    fastas = q_runner(n_procs, cluster, worker1, conf, tdir)
    print 'past fastas process'
    #!! RepeatMasker with multicore options (so no q_runner or runner action)
    repeatmasker(n_procs, fastas)
    fstas = [f[1] for f in fastas]
    print 'fstas =', fstas
    #pdb.set_trace()
    # convert the .mask files to 2 bit
    # !!! decided not to do this because faToTwoBit acts weird with subprocess
    #fstas = ['/tmp/msatcommander-454-fACVUg/clean-moss_2ee-wytONi.fa']
    #pdb.set_trace()
    twobits = q_runner(n_procs, fstas, worker2)
    print 'twobits=', twobits
    print 'past twobit process'
    #pdb.set_trace()
    # TODO:  cleanup tmp a little to conserve space
    m_data = q_runner(n_procs, cluster, worker1, conf, tdir, True)
    #pdb.set_trace()
    msats = [m[1] for m in m_data]
    print 'msat =', msats
    # we need to zip the fa/2bits and the msats together, but they could be out
    # of order.  so:
    file_stack = ()
    for c in cluster:
        for t in fstas:
            if c in t:
                tb = t
        for m in msats:
            if c in m:
                mb = m
        file_stack = file_stack + ((tb, mb,),)
    #pdb.set_trace()
    #print 'waiting 5s'
    #time.sleep(5)
    print 'file_stack= ', file_stack
    print 'starting blat'
    blat = q_runner(n_procs, file_stack, worker3, tdir)
    pdb.set_trace()
    for m in blat:
        chunks = list(chunker(m[1]))
        parser = q_runner(n_procs, chunks, worker4, conf)
            
    
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'
    #pdb.set_trace()


if __name__ == '__main__':
    main()

