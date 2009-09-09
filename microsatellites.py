#!/usr/bin/env python
# encoding: utf-8
"""
main2.py

Created by Brant Faircloth on 2009-09-07.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import re, pdb, time, MySQLdb, ConfigParser, multiprocessing, cPickle, msat, new, numpy
import MySQLdb.cursors


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

def worker(id, record, motifs):
    # find msats
    record = msatSearch(record, motifs)
    # mask msat sequence
    record = softmask(record)
    # open connection to microsatellites table
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
    db="454_msatcommander")
    cur = conn.cursor()
    # add new data for mask table (which = microsatellite repeats)
    for match in record.matches:
        for repeat in record.matches[match]:
            motif_count = (repeat[0][1] - repeat[0][0])/len(match)
            #pdb.set_trace()
            cur.execute('''INSERT INTO mask (name, motif, start, end, 
            preceding, following, motif_count) VALUES (%s,%s,%s,%s,%s,%s,%s) 
            ''', (record.id, match, repeat[0][0], repeat[0][1], repeat[1], 
            repeat[2], motif_count))
    # update blob in 'main' table
    if record.matches:
        record_pickle = cPickle.dumps(record,1)
        # replace sequence with masked version in seq_trimmed - we can always
        # call UPPER() on it for non-masked version
        cur.execute('''UPDATE sequence_test SET seq_trimmed = %s, record = %s, msat = %s WHERE id = %s''', (record.seq.masked, record_pickle, True, id))
    else:
        cur.execute('''UPDATE sequence_test SET msat = %s WHERE id = %s''', (False, id))
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

def createMaskTable(cur):
    try:
        cur.execute('''DROP TABLE mask''')
    except:
        pass
    cur.execute('''CREATE TABLE mask (id INT UNSIGNED NOT NULL 
            AUTO_INCREMENT, name VARCHAR(100), motif VARCHAR(8), start 
            SMALLINT UNSIGNED, end SMALLINT UNSIGNED, preceding SMALLINT 
            UNSIGNED, following SMALLINT UNSIGNED, motif_count SMALLINT 
            UNSIGNED, PRIMARY KEY (id), INDEX mask_name (name))''')

def updateSequenceTable(cur):
    try:
        cur.execute('''ALTER TABLE sequence_test ADD COLUMN msat BOOLEAN''')
    except:
        pass

def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    # get sequence from dbase.  Use the server-side cursor to keep the load
    # down from large result sets
    #conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
    #db="454_msatcommander", cursorclass=MySQLdb.cursors.SSCursor)
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
    db="454_msatcommander")
    cur = conn.cursor()
    createMaskTable(cur)
    cur.execute('''SELECT id, record FROM sequence_test WHERE n_count <= 2 AND
    trimmed_len > 10''')
    # pdb.set_trace()
    # setup motifs
    motifs = motifCollection(min_length = [10,6,4,4,4,4], scan_type = "2+", \
    perfect = True)
    # add microsatellite method to sequence
    if conf.getboolean('Multiprocessing', 'MULTIPROCESSING'):
        # get num processors
        n_procs = conf.get('Multiprocessing','processors')
        if n_procs == 'Auto':
            n_procs = multiprocessing.cpu_count() - 1
        else:
            n_procs = int(n_procs)
        print 'Multiprocessing.  Number of processors = ', n_procs
        # to test with fewer sequences
        #count = 0
        threads = []
        # access the data on sequence by sequence basis to avoid reading the 
        # entire table contents into memory
        row = cur.fetchone()
        while row is not None:
        #index = 0
        #while index < 5000:
            if len(threads) < n_procs:
                #row = cur.fetchone()
                # convert BLOB back to sequence record
                id, record = row[0], cPickle.loads(row[1])
                p = multiprocessing.Process(target=worker, args=(id, record, motifs))
                p.start()
                threads.append(p)
                row = cur.fetchone()
                #index += 1
            else:
                for t in threads:
                    if not t.is_alive():
                        threads.remove(t)
    else:
        print 'Not using multiprocessing'
        # access the data on sequence by sequence basis to avoid 
        # reading the entire table contents into memory
        row = cur.fetchone()
        while row is not None:
            # convert BLOB back to sequence record
            id, record = row[0], cPickle.loads(row[1])
            pdb.set_trace()
            worker(record, motifs)
            row = cur.fetchone()
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'


if __name__ == '__main__':
    main()

