#!/usr/bin/env python
# encoding: utf-8
"""
main2.py

Created by Brant Faircloth on 2009-09-07.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os
import sys
import re
import pdb
import time
import MySQLdb
import ConfigParser
import multiprocessing
import cPickle
#import new
import numpy
import progress
import optparse
import operator
import MySQLdb.cursors
import msat


def softmask(record):
    sequence    = numpy.array(list(str(record.seq)))
    lower       = numpy.array(list(str(record.seq).lower()))
    mask        = numpy.ma.make_mask_none((len(sequence),))
    for match in record.matches.values():
        for repeat in match:
            mask[repeat[0][0]:repeat[0][1]] = True
    # we're using a masking array to soft-mask the sequences versus iterating
    # over the bases
    numpy.putmask(sequence, mask, lower)
    record.seq.masked = sequence.tostring()
    return record


def microsatellite(record, msat):
    '''Generalized microsatellite search function'''
    #pdb.set_trace()
    for repeat in range(len(msat.compiled)):
        temp_match = ()
        for m in msat.compiled[repeat].finditer(str(record.seq)):
            temp_match += ((m.span(),m.span()[0],len(record.seq)-m.span()[1]),)
        if temp_match:
            record.matches[msat.motif[repeat]] = temp_match

def msatSearch(record, motifs):
    # add matches attribute to object
    record.matches = {}
    # add search method to object
    #record.microsatellite = \
    #new.instancemethod(microsatellite, record, \
    #record.__class__)
    for search in motifs:
        microsatellite(record, search)
    return record

def combineLoci(record, min_distance):
    '''combined adjacent loci - this is somewhat cumbersome due to the
    format of the matches returned from msat (a dict with keys = motif).
    Essentially, we are running a pairwise comparison across all motifs
    located to determine which are within a predetermined distance from 
    one another'''
    record.combined = ()
    if record.matches:
        temp_combined = []
        reorder = ()
        # turn our dict into something more useful for this purpose
        for motif in record.matches:
            for pos,val in enumerate(record.matches[motif]):
                reorder += ((motif, pos, val[0][0], val[0][1], val[1], val[2]),)
        # sort it
        reorder = sorted(reorder, key=operator.itemgetter(2))
        # combine adjacent loci at < min_distance
        for k,v in enumerate(reorder):
            included = False
            if not temp_combined:
                temp_combined.append([v])
            else:
                for gp, g in enumerate(temp_combined):
                    for elem in g:
                        if v[2] - elem[3] <= min_distance:
                            temp_combined[gp].append(v)
                            included = True
                            break
                # ensure we add those that do not combine
                if not included:
                    temp_combined.append([v])
        for group in temp_combined:
            motifs = []
            if len(group) > 1:
                gs = group[0][2]
                ge = group[-1][3]
                gp = group[0][4]
                gf = group[-1][5]
            else:
                gs, ge = group[0][2], group[0][3]
                gp, gf = group[0][4], group[0][5]
            name = ''
            member_count = 0
            for pos,member in enumerate(group):
                if pos + 1 < len(group):
                    dist = group[pos + 1][3] - group[pos][3]
                    if dist > 1:
                        spacer = '...'
                    else:
                        spacer = ''
                else:
                    spacer = ''
                length = (member[3]-member[2])/len(member[0])
                name += '%s(%s)%s' % (member[0], length, spacer)
                motifs.append([member[0],length])
                member_count += 1
            record.combined += (((gs, ge), gp, gf, member_count, motifs, name),)
    return record

def worker(id, record, motifs, db, have_sequence_table, combine_loci, combine_loci_dist):
    # find msats
    record = msatSearch(record, motifs)
    # only mask msat sequence if we started with a sequence table
    if have_sequence_table:
        record = softmask(record)
    # combined loci
    if combine_loci:
        record = combineLoci(record, combine_loci_dist)
        #pdb.set_trace()
    # open connection to microsatellites table
    conn = MySQLdb.connect(user     = db[0], 
                           passwd   = db[1], 
                           db       = db[2]
                           )
    cur = conn.cursor()
    # add new data for mask table (which = microsatellite repeats)
    # setup our own auto-incrementing index, so we can use value
    # later without multiprocessing causing us a problem
    msat_id = 0
    for match in record.matches:
        for repeat in record.matches[match]:
            motif_count = (repeat[0][1] - repeat[0][0])/len(match)
            #pdb.set_trace()
            cur.execute('''INSERT INTO mask (sequence_id, id, motif, start, end, 
            preceding, following, motif_count) VALUES (%s, %s, %s, %s, %s, %s, %s, %s) 
            ''', (id, msat_id, match, repeat[0][0], repeat[0][1], repeat[1], 
            repeat[2], motif_count))
            msat_id += 1
    if combine_loci:
        # setup our own auto-incrementing index, so we can use value
        # later without multiprocessing causing us a problem
        combined_id = 0
        for combined in record.combined:                
            cur.execute('''INSERT INTO combined (sequence_id, id, motif, start, end, 
            preceding, following, members) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ''', (id, combined_id, combined[-1], combined[0][0], combined[0][1], combined[1], 
            combined[2], combined[3]))
            for m in combined[4]:
                cur.execute('''INSERT INTO combined_components 
                    (sequence_id, combined_id, motif, length) VALUES 
                    (%s, %s, %s, %s)''', (id, combined_id, m[0], m[1]))
            combined_id += 1
    # update blob in 'main' table
    if have_sequence_table:
        if record.matches:
            record_pickle = cPickle.dumps(record,1)
            # replace sequence with masked version in seq_trimmed - we can always
            # call UPPER() on it for non-masked version
            cur.execute('''UPDATE sequence SET seq_trimmed = %s, 
                record = %s, msat = %s WHERE id = %s''', 
                (record.seq.masked, record_pickle, True, id))
        else:
            cur.execute('''UPDATE sequence SET msat = %s WHERE id = %s''', 
                (False, id))
    cur.close()
    conn.commit()
    conn.close()
    
def createMotifInstances(motif, min_length, perfect):
    return msat.seqsearch.MicrosatelliteMotif(motif, min_length, perfect)

def motifCollection(**kwargs):
    possible_motifs = (msat.motif.mononucleotide, msat.motif.dinucleotide, \
    msat.motif.trinucleotide, msat.motif.tetranucleotide,
    msat.motif.pentanucleotide, msat.motif.hexanucleotide)
    # add optional lengths
    possible_motifs = zip(kwargs['min_length'], possible_motifs)
    collection = ()
    if kwargs['scan_type'] == 'all':
        for m in possible_motifs:
            pdb.set_trace()
            collection += (createMotifInstances(m[1], m[0], \
            kwargs['perfect']),)
    elif '+' in kwargs['scan_type']:
        # subtracting 1 so that we get >= options.scan_type
        scan = int(kwargs['scan_type'][0]) - 1
        for m in possible_motifs[scan:]:
            collection += (createMotifInstances(m[1], m[0], \
            kwargs['perfect']),)
    elif '-' in kwargs['scan_type']:
        scan_start = int(kwargs['scan_type'][0]) - 1
        scan_stop = int(kwargs['scan_type'][2])
        for m in possible_motifs[scan_start:scan_stop]:
            collection += (createMotifInstances(m[1], m[0], \
            kwargs['perfect']),)
    else:
        # no iteration here because tuple != nested
        scan = int(kwargs['scan_type'][0]) - 1
        collection += (createMotifInstances(possible_motifs[scan][1], \
        possible_motifs[scan][0], kwargs['perfect']),)
    return collection


def createSequenceTable(cur):
    """create a quasi-sequence table to hold identifying information for
    each fasta or 2bit read - used only where SequenceTable = False"""
    try:
        cur.execute('''DROP TABLE sequence''')
    except:
        pass
    cur.execute(''' CREATE TABLE sequence (
    id int(10) unsigned NOT NULL,
    name varchar(100) DEFAULT NULL,
    PRIMARY KEY (id),
    KEY name (name)
    ) ENGINE=InnoDB''')
    

def createMaskTableWithForeign(cur):
    try:
        cur.execute('''DROP TABLE mask''')
    except:
        pass
    # TODO:  Switch index to reference sequence.id versus autoincrement value and create an index on it
    cur.execute('''CREATE TABLE mask (
        sequence_id INT(10) UNSIGNED NOT NULL, 
        id int(10) UNSIGNED NOT NULL, 
        motif VARCHAR(8), 
        start BIGINT UNSIGNED, 
        end BIGINT UNSIGNED, 
        preceding BIGINT UNSIGNED, 
        following BIGINT UNSIGNED, 
        motif_count MEDIUMINT UNSIGNED, 
        PRIMARY KEY(sequence_id, id),
        FOREIGN KEY (sequence_id) REFERENCES sequence (id)) ENGINE=InnoDB''')

def createCombinedLociWithForeign(cur):
    try:
        cur.execute('''DROP TABLE combined''')
    except:
        pass
    # TODO:  Switch index to reference sequence.id versus autoincrement value and create an index on it
    cur.execute('''CREATE TABLE combined (
        sequence_id INT(10) UNSIGNED NOT NULL, 
        id int(10) UNSIGNED NOT NULL,
        motif TEXT, 
        start BIGINT UNSIGNED, 
        end BIGINT UNSIGNED,
        preceding BIGINT UNSIGNED,
        following BIGINT UNSIGNED,
        members MEDIUMINT UNSIGNED,
        PRIMARY KEY(sequence_id, id),
        INDEX(members),
        INDEX(id),
        FOREIGN KEY (sequence_id) REFERENCES sequence (id)) ENGINE=InnoDB''')
    createCombinedComponentsTable(cur)

def createCombinedComponentsTable(cur):
    try:
        cur.execute('''DROP TABLE combined_components''')
    except:
        pass
    cur.execute('''CREATE TABLE combined_components (
        sequence_id INT(10) UNSIGNED NOT NULL,
        combined_id INT(10) UNSIGNED NOT NULL,
        motif TEXT,
        length MEDIUMINT UNSIGNED NOT NULL,
        INDEX(combined_id),
        FOREIGN KEY (sequence_id, combined_id) REFERENCES combined (sequence_id, id)
        ) ENGINE=InnoDB''')


def updateSequenceTable(cur):
    # TODO:  This is dumb.  Just add the column in the linkers.py
    try:
        cur.execute('''ALTER TABLE sequence ADD COLUMN msat BOOLEAN''')
    except MySQLdb._mysql.OperationalError, e:
        if e[0] == 1060:
            cur.execute('''UPDATE sequence SET msat = NULL''')
            print 'Zeroed-out sequence.msat column'

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
    # Copyright (c) 2010 Brant C. Faircloth                      #
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

class Sequence():
    """docstring for Sequence"""
    def __init__(self, engine='mysql', **kwargs):
        self.engine = engine
        if self.engine == 'mysql':
            self.create_mysql_iterator(**kwargs)
        elif self.engine == 'biopython' and kwargs['data_type'] == 'fasta':
            self.create_biopython_iterator(**kwargs)
        elif self.engine == 'pyfasta' and kwargs['data_type'] == 'fasta':
            self.create_pyfasta_iterator(**kwargs)
        elif self.engine == 'twobit' and kwargs['data_type'] == 'twobit':
            self.create_twobit_iterator(**kwargs)
    
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
        import bx.seq.twobit
        self.fasta = bx.seq.twobit.TwoBitFile(file(kwargs['input']))
        self.readcount = self.fasta.seq_count
        self.db_values = zip(range(self.fasta.seq_count), sorted(self.fasta.keys()))
        self.read = iter(self.db_values)
    
    def create_pyfasta_iterator(self, **kwargs):
        from pyfasta import Fasta
        print "Generating PyFasta sequence index.  This may take a moment...."
        self.fasta = Fasta(kwargs['input'])
        self.readcount = len(self.fasta)
        self.db_values = zip(range(len(self.fasta)), sorted(self.fasta.keys()))
        self.read = iter(self.db_values)

class SequenceWrapper():
    """this class wraps pyfasta objects to make them appear similar to 
    biopython sequence records"""
    def __init__(self, iden, sequence):
        self.id = iden
        self.seq = sequence

def drop_old_tables(cur, have_sequence_table):
    for t in ['combined_components', 'combined', 'mask']:
        query = 'DROP TABLE IF EXISTS {0}'.format(t)
        cur.execute(query)
    if not have_sequence_table:
        cur.execute('''DROP table IF EXISTS sequence''')

def main():
    start_time = time.time()
    options, arg = interface()
    motd()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    # =============================================
    # = Setup additional configuration parameters =
    # =============================================
    db = (conf.get('Database','USER'), 
          conf.get('Database','PASSWORD'), 
          conf.get('Database','DATABASE')
          )
    have_sequence_table = conf.getboolean('MicrosatelliteParameters', 'HaveSequenceTable')
    fasta_engine = conf.get('MicrosatelliteParameters', 'FastaEngine').lower()
    try:
        data_type = conf.get('Input', 'Type')
    except:
        data_type = 'fasta'
    combine_loci = conf.getboolean('MicrosatelliteParameters', 'CombineLoci')
    combine_loci_dist = conf.getint('MicrosatelliteParameters', 'CombineLociDist')
    m_processing = conf.getboolean('Multiprocessing', 'MULTIPROCESSING')
    get_num_procs = conf.get('Multiprocessing','processors')
    
    conn = MySQLdb.connect(user     = db[0],
                           passwd   = db[1],
                           db       = db[2]
                           )
    cur = conn.cursor()
    # Drop old tables
    drop_old_tables(cur, have_sequence_table)
    if have_sequence_table:
        data = Sequence(engine = 'mysql', cursor = cur)
    else:
        # create a quasi-sequence table
        createSequenceTable(cur)
        # get out data
        data = Sequence(engine = fasta_engine, 
            input = conf.get('Input','sequence'), 
            data_type = data_type
            )
        cur.executemany('''INSERT INTO sequence (id, name) VALUES (%s, %s)''', data.db_values)
    # create our msat table
    createMaskTableWithForeign(cur)
    # create the combined msat table
    if combine_loci:
        createCombinedLociWithForeign(cur)
    conn.commit()
    scan_type = conf.get('MicrosatelliteParameters', 'ScanType')
    motifs = motifCollection(min_length = [10,6,4,4,4,4], scan_type = scan_type, \
                perfect = True)
    if m_processing:
        # get num processors
        n_procs = get_num_procs
        if n_procs == 'Auto':
            n_procs = multiprocessing.cpu_count() - 2
        else:
            n_procs = int(n_procs)
        print 'Multiprocessing.  Number of processors = %s\n' % n_procs
        # to test with fewer sequences
        #count = 0
        threads = []
        # access the data on sequence by sequence basis to avoid reading the 
        # entire table contents into memory        
        pb = progress.bar(0,data.readcount,60)
        pb_inc = 0
        try:
            while data:
                if len(threads) < n_procs:
                    # convert BLOB back to sequence record
                    if data.engine == 'mysql':
                        # convert BLOB back to sequence record
                        record = data.read.next()
                        iden = record[0]
                        record = cPickle.loads(record[1])
                    elif data.engine == 'pyfasta' or data.engine == 'biopython':
                        iden, chromo = data.read.next()
                        record = SequenceWrapper(iden, data.fasta[chromo])
                    elif data.engine == 'twobit':
                        iden, chromo = data.read.next()
                        record = SequenceWrapper(iden, data.fasta[chromo][:])
                    p = multiprocessing.Process(target=worker, args=(
                                    iden,
                                    record,
                                    motifs,
                                    db,
                                    have_sequence_table,
                                    combine_loci, 
                                    combine_loci_dist)
                                    )
                    p.start()
                    threads.append(p)
                    if (pb_inc+1)%1000 == 0:
                        pb.__call__(pb_inc+1)
                    elif pb_inc + 1 == data.readcount:
                        pb.__call__(pb_inc+1)
                    pb_inc += 1
                else:
                    for t in threads:
                        if not t.is_alive():
                            threads.remove(t)
        except StopIteration:
            pass
    else:
        print 'Not using multiprocessing\n'
        # access the data on sequence by sequence basis to avoid 
        # reading the entire table contents into memory
        pb = progress.bar(0,data.readcount,60)
        pb_inc = 0
        #pdb.set_trace()
        try:
            #pdb.set_trace()
            while data:
                if data.engine == 'mysql':
                    # convert BLOB back to sequence record
                    record = data.read.next()
                    iden = record[0]
                    record = cPickle.loads(record[1])
                elif data.engine == 'pyfasta':
                    iden, chromo = data.read.next()
                    record = SequenceWrapper(iden, data.fasta[chromo])
                elif data.engine == 'biopython':
                    iden, chromo = data.read.next()
                    record = data.fasta[chromo]
                elif data.engine == 'twobit':
                    iden, chromo = data.read.next()
                    record = SequenceWrapper(iden, data.fasta[chromo][:])
                #pdb.set_trace()
                worker(iden, record, motifs, db, have_sequence_table,
                        combine_loci, combine_loci_dist)
                row = cur.fetchone()
                if (pb_inc+1)%1000 == 0:
                    pb.__call__(pb_inc+1)
                elif pb_inc + 1 == data.readcount:
                    pb.__call__(pb_inc+1)
                pb_inc += 1
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

