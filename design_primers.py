#!/usr/bin/env python
# encoding: utf-8

"""
design_primers.py

Created by Brant Faircloth on 10 May 2010 08:47 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import os
import sys
import time
import string
import MySQLdb
import cPickle
import optparse
import progress
import ConfigParser
import multiprocessing
import p3wrapr


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', dest = 'conf', action='store', 
type='string', default = None, help='The path to the configuration file.', 
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

def create_primers_table(cur):
    try:
        cur.execute('''DROP TABLE primers''')
    except:
        pass
    cur.execute('''CREATE TABLE primers (
        sequence_id int(10) UNSIGNED NOT NULL,
        id int(10) UNSIGNED NOT NULL,
        primer int(4) UNSIGNED NOT NULL,
        left_p varchar(100) NOT NULL,
        left_sequence varchar(100) NOT NULL,
        left_tm float,
        left_gc float,
        left_self_end float,
        left_self_any float,
        left_hairpin float,
        left_end_stability float,
        left_penalty float,
        right_p varchar(100) NOT NULL,
        right_sequence varchar(100) NOT NULL,
        right_tm float,
        right_gc float,
        right_self_end float,
        right_self_any float,
        right_hairpin float,
        right_end_stability float,
        right_penalty float,
        pair_product_size float,
        pair_compl_end float,
        pair_compl_any float,
        pair_penalty float,
        INDEX seq_idx (sequence_id, id),
        FOREIGN KEY(sequence_id, id) REFERENCES combined(sequence_id, id)
        ) ENGINE=InnoDB''')

def create_tagged_primers_table(cur):
    try:
        cur.execute('''DROP TABLE tagged_primers''')
    except:
        pass
    cur.execute('''CREATE TABLE tagged_primers (
        sequence_id int(10) UNSIGNED NOT NULL,
        id int(10) UNSIGNED NOT NULL,
        primer int(4) UNSIGNED NOT NULL,
        best int(4) UNSIGNED NOT NULL,
        tag varchar(100),
        tagged varchar(100),
        tag_seq varchar(100),
        common varchar(100),
        pigtail_tagged varchar(100),
        pigtail_tag_seq varchar(100),
        pigtail_common varchar(100),
        left_p varchar(100),
        left_sequence varchar(100),
        left_self_end float,
        left_self_any float,
        left_hairpin float,
        left_penalty float,
        right_p varchar(100),
        right_sequence varchar(100),
        right_self_end float,
        right_self_any float,
        right_hairpin float,
        right_penalty float,
        pair_product_size float,
        pair_compl_end float,
        pair_compl_any float,
        pair_penalty float,
        INDEX seq_idx (sequence_id, id),
        FOREIGN KEY(sequence_id, id) REFERENCES combined(sequence_id, id)
        ) ENGINE=InnoDB''')


def insert_primers(cur, seq_iden, iden, absolute_start, primer3):
    for i,p in primer3.primers.iteritems():
        if i != 'metadata':
            # create a copy of the dict, to which we add the
            # FOREIGN KEY reference
            td = p.copy()
            td['SEQUENCE_ID'] = seq_iden
            td['ID'] = iden
            td['PRIMER'] = i
            for pos in ['PRIMER_LEFT', 'PRIMER_RIGHT']:
                start, end = td[pos].split(',')
                td[pos] = "{0},{1}".format(absolute_start + int(pos), end)
            cur.execute('''INSERT INTO primers VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
                )''', (
                td['SEQUENCE_ID'],
                td['ID'],
                td['PRIMER'],
                td['PRIMER_LEFT'],
                td['PRIMER_LEFT_SEQUENCE'],
                td['PRIMER_LEFT_TM'],
                td['PRIMER_LEFT_GC_PERCENT'],
                td['PRIMER_LEFT_SELF_END_TH'],
                td['PRIMER_LEFT_SELF_ANY_TH'],
                td['PRIMER_LEFT_HAIRPIN_TH'],
                td['PRIMER_LEFT_END_STABILITY'],
                td['PRIMER_LEFT_PENALTY'],
                td['PRIMER_RIGHT'],
                td['PRIMER_RIGHT_SEQUENCE'],
                td['PRIMER_RIGHT_TM'],
                td['PRIMER_RIGHT_GC_PERCENT'],
                td['PRIMER_RIGHT_SELF_END_TH'],
                td['PRIMER_RIGHT_SELF_ANY_TH'],
                td['PRIMER_RIGHT_HAIRPIN_TH'],
                td['PRIMER_RIGHT_END_STABILITY'],
                td['PRIMER_RIGHT_PENALTY'],
                td['PRIMER_PAIR_PRODUCT_SIZE'],
                td['PRIMER_PAIR_COMPL_END_TH'],
                td['PRIMER_PAIR_COMPL_ANY_TH'],
                td['PRIMER_PAIR_PENALTY']))

def insert_tagged_primers(cur, conn, seq_iden, iden, primer3):
    if primer3.tagged_good:
        for i,p in primer3.tagged_good.iteritems():
            if i != 'metadata':
                # create a copy of the dict, to which we add the
                # FOREIGN KEY reference
                td = p.copy()
                td['SEQUENCE_ID'] = seq_iden
                td['ID'] = iden
                td['BEST']      = 0
                td['PRIMER'], td['TAG'], td['TAGGED'] = i.split('_')
                cur.execute('''INSERT INTO tagged_primers VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)''', (
                    td['SEQUENCE_ID'],
                    td['ID'],
                    td['PRIMER'],
                    td['BEST'],
                    td['TAG'],
                    td['PRIMER_TAGGED'],
                    td['PRIMER_TAG'],
                    td['PRIMER_TAG_COMMON_BASES'],
                    td['PRIMER_PIGTAILED'],
                    td['PRIMER_PIGTAIL_TAG'],
                    td['PRIMER_PIGTAIL_TAG_COMMON_BASES'],
                    td['PRIMER_LEFT'],
                    td['PRIMER_LEFT_SEQUENCE'],
                    td['PRIMER_LEFT_SELF_END_TH'],
                    td['PRIMER_LEFT_SELF_ANY_TH'],
                    td['PRIMER_LEFT_HAIRPIN_TH'],
                    td['PRIMER_LEFT_PENALTY'],
                    td['PRIMER_RIGHT'],
                    td['PRIMER_RIGHT_SEQUENCE'],
                    td['PRIMER_RIGHT_SELF_END_TH'],
                    td['PRIMER_RIGHT_SELF_ANY_TH'],
                    td['PRIMER_RIGHT_HAIRPIN_TH'],
                    td['PRIMER_RIGHT_PENALTY'],
                    td['PRIMER_PAIR_PRODUCT_SIZE'],
                    td['PRIMER_PAIR_COMPL_END_TH'],
                    td['PRIMER_PAIR_COMPL_ANY_TH'],
                    td['PRIMER_PAIR_PENALTY'])
                    )
                #pdb.set_trace()
        conn.commit()
        if primer3.tagged_best:
            best = primer3.tagged_best.keys()[0].split('_')
            # we're using real side names now
            if best[2]=='r':
                side = 'RIGHT'
            else:
                side = 'LEFT'
            #pdb.set_trace()
            cur.execute('''UPDATE tagged_primers 
                SET best = 1 WHERE 
                sequence_id = %s
                AND id = %s
                AND primer = %s
                AND tag = %s 
                AND tagged = %s''',
                (seq_iden, iden, best[0], best[1], side))
        conn.commit()


class Sequence():
    """docstring for Sequence"""
    def __init__(self, engine='mysql', function = 'iterator', **kwargs):
        self.engine = engine
        if self.engine == 'mysql' and function == 'iterator':
            self.create_mysql_iterator(**kwargs)
        elif self.engine == 'biopython' and kwargs['data_type'] == 'fasta' and function == 'iterator':
            self.create_biopython_iterator(**kwargs)
        elif self.engine == 'pyfasta' and function == 'iterator':
            self.create_pyfasta_iterator(**kwargs)
        elif self.engine =='pyfasta' and function == 'getter':
            self.get_pyfasta_reads(**kwargs)

    def create_mysql_iterator(self, **kwargs):
        cur = kwargs['cursor']
        query = '''SELECT id, record FROM sequence WHERE n_count <= 2 AND 
                    trimmed_len > 40'''
        cur.execute(query)
        self.readcount = cur.rowcount
        self.read = iter(cur.fetchall())

    def create_biopython_iterator(self, **kwargs):
        from Bio import SeqIO
        print "Generating BioPython sequence index.  This may take a moment...."
        self.fasta = SeqIO.index(kwargs['input'], kwargs['data_type'])
        self.readcount = len(self.fasta)
        self.db_values = zip(range(len(self.fasta)), sorted(self.fasta.keys()))
        self.read = iter(self.db_values)

    def create_twobit_iterator(self, **kwargs):
        pass

    def create_pyfasta_iterator(self, **kwargs):
        from pyfasta import Fasta
        print "Generating PyFasta sequence index.  This may take a moment...."
        self.fasta = Fasta(kwargs['input'])
        self.readcount = len(self.fasta)
        self.db_values = zip(range(len(self.fasta)), sorted(self.fasta.keys()))
        self.read = iter(self.db_values)

    def get_pyfasta_reads(self, **kwargs):
        from pyfasta import Fasta
        self.fasta = Fasta(kwargs['input'])
        self.readcount = len(self.fasta)



def get_sequence_for_chunk(data, chunk, flank = 300):
    # maketrans is _way_ faster than re.sub
    trans = string.maketrans('acgt', 'NNNN')
    container = []
    for row in chunk:
        # flank out 300 bp in each direction - do this the "hard" way
        # so that we can deal with any faults where the snip on either side
        # extends beyond the length of the physical sequence.  In the current
        # case we should fail gracefully, and still attempt to design primers
        # with what we've got
        chromo = data.fasta[row[0]]
        start = (row[3] - flank)
        end = (row[4] + flank)
        if start < 0 and end > len(chromo):
            preceding = chromo[0:row[3]]
            middle = chromo[row[3]:row[4]]
            following = chromo[row[4]:len(chromo)]
            sequence.absolute_start = 0
        
        elif start < 0:
            preceding = chromo[0:row[3]]
            middle = chromo[row[3]:row[4]]
            following = chromo[row[4]:end]
            sequence.absolute_start = 0
        
        elif end > len(data.fasta[row[0]]):
            preceding = chromo[start:row[3]]
            middle = chromo[row[3]:row[4]]
            following = chromo[row[4]:len(chromo)]
            sequence.absolute_start = start
        
        else:
            preceding = chromo[start:row[3]]
            middle = chromo[row[3]:row[4]]
            following = following = chromo[row[4]:end]
            sequence.absolute_start = start
            
        temp_seq = preceding + middle + following
        # convert masked bases to N
        mask_seq = string.translate(temp_seq, trans)
        iden = '{0}_{1}_{2}'.format(row[0],row[1],row[2])
        sequence = SequenceWrapper(iden, mask_seq)
        sequence.chromo = row[0]
        sequence.sequence_id = row[1]
        sequence.id = row[2]
        # get relative start and end positions
        sequence.repeat_relative_start = len(preceding)
        sequence.repeat_relative_end = len(preceding) + len(middle)
        # add the sequences to the container
        container.append(sequence)
    return container

class SequenceWrapper():
    """this class wraps pyfasta objects to make them appear similar to 
    biopython sequence records"""
    def __init__(self, iden, sequence):
        self.id = iden
        self.seq = sequence

def worker(db, container, settings, tag_settings):
    """docstring for worker"""
    # open our own connection to the db
    conn = MySQLdb.connect(
        user    = db[0],
        passwd  = db[1],
        db      = db[2]
        )
    cur = conn.cursor()
    for sequence in container:
        target = '{0},{1}'.format(sequence.repeat_relative_start, 
            sequence.repeat_relative_end - sequence.repeat_relative_start)
        primer3 = p3wrapr.Primers()
        primer3.pick(settings, sequence = str(sequence.seq), target=target, name = 'primers')
        #primer3.tag(tag_settings, M13F='GTAAAACGACGGCCAG')
        primer3.tag(tag_settings, CAG = 'CAGTCGGGCGTCATCA', M13R = 'GGAAACAGCTATGACCAT')
        if primer3.tagged_good:
            primer3.pigtail(tag_settings, 'GTTT')
        insert_primers(cur, sequence.sequence_id, sequence.id, sequence.absolute_start, primer3)
        insert_tagged_primers(cur, conn, sequence.sequence_id, sequence.id, primer3)
    # make sure we commit our results before closing
    conn.commit()
    # close our db connection
    cur.close()
    conn.close()
    


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    options, arg    = interface()
    conf            = ConfigParser.ConfigParser()
    conf.read(options.conf)
    have_sequence_table = conf.getboolean('MicrosatelliteParameters', 'HaveSequenceTable')
    db = (conf.get('Database','USER'), 
          conf.get('Database','PASSWORD'), 
          conf.get('Database','DATABASE')
          )
    m_processing = conf.getboolean('Multiprocessing', 'MULTIPROCESSING')
    num_procs = conf.get('Multiprocessing','processors')
    # build our configuration
    conn = MySQLdb.connect(
        user    = db[0],
        passwd  = db[1],
        db      = db[2]
        )
    cur = conn.cursor()
    # get groups of 1000 sequences, package them, send them to muliprocessing
    # for primer design
    db_chunk = []
    # get count of total combined msats in the table
    cur.execute('''SELECT sequence.name, combined.sequence_id, combined.id, 
            combined.start, combined.end FROM sequence, combined 
            where combined.sequence_id = sequence.id''')
    sequence_ids = cur.fetchall()
    # split our msats into chunks or blocks of 1000 reads, one of which we
    # will sequentially pass to each processing core across the set of reads
    if len(sequence_ids) > 1000:
        for span in range(0,len(sequence_ids), 1000):
            db_chunk.append(sequence_ids[span:span + 1000])
    db_chunk_len = float(len(db_chunk))
    db_chunk = iter(db_chunk)
    data = Sequence(input = conf.get('Input','sequence'), engine = 'pyfasta', function = 'getter')
    # we only need a single instance of our primer design settings
    settings = p3wrapr.primer.Settings()
    settings.basic()
    # ensure that we get pretty different primers for each set of X designed
    # => offset by 10 bp.  this makes primer3 somewhat slower.
    settings.params['PRIMER_MIN_THREE_PRIME_DISTANCE'] = 10
    # override the default location for the mispriming library
    settings.params['PRIMER_MISPRIMING_LIBRARY'] = '/Users/bcf/Bin/misprime_lib_weight'
    # we only need a single instance of our tag settings
    tag_settings = p3wrapr.primer.Settings()
    tag_settings.reduced(PRIMER_PICK_ANYWAY=1)
    # create primers and tagged primers table
    create_primers_table(cur)
    create_tagged_primers_table(cur)
    if m_processing:
        if num_procs == 'Auto':
            num_procs = multiprocessing.cpu_count() - 2
        else:
            num_procs = int(num_procs)
        print 'Multiprocessing.  Number of processors = %s\n' % num_procs
        threads = []
        # access the data on sequence by sequence basis to avoid reading the 
        # entire table contents into memory        
        pb = progress.bar(0,db_chunk_len,60)
        pb_inc = 0
        try:
            while db_chunk:
                if len(threads) < num_procs:
                    container = get_sequence_for_chunk(data, db_chunk.next())
                    p = multiprocessing.Process(target=worker, args=(db, container, settings, tag_settings))
                    p.start()
                    threads.append(p)
                    if (pb_inc+1)%1 == 0:
                        pb.__call__(pb_inc+1)
                    elif pb_inc + 1 == db_chunk_len:
                        pb.__call__(pb_inc+1)
                    pb_inc += 1
                else:
                    for t in threads:
                        if not t.is_alive():
                            threads.remove(t)
        except StopIteration:
            pass
    # for single processing
    else:
        try:
            while db_chunk:
                if have_sequence_table:
                    pass
                else:
                    container = get_sequence_for_chunk(data, db_chunk.next())
                    worker(db, container, settings, tag_settings)
        except StopIteration:
            pass
    print '\n'
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()