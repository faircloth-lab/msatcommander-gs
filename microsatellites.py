#!/usr/bin/env python
# encoding: utf-8
"""
main2.py

Created by Brant Faircloth on 2009-09-07.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os, sys, re, pdb, time, MySQLdb, ConfigParser, multiprocessing, cPickle, msat, new, numpy, progress, optparse, operator
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

def combineLoci(record, min_distance):
    '''combined adjacent loci - this is somewhat cumbersome due to the
    format of the matches returned from msat (a dict with keys = motif).
    Essentially, we are running a pairwise comparison across all motifs
    located to determine which are within a predetermined distance from 
    one another'''
    record.combined = {}
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
        for i in reorder:
            included = False
            if not temp_combined:
                temp_combined.append([i])
            else:
                for gp, g in enumerate(temp_combined):
                    for elem in g:
                        if i[2] - elem[3] <= min_distance:
                            temp_combined[gp].append(i)
                            included = True
                            break
                if not included:
                    temp_combined.append([i])
        #if len(temp_combined[0]) > 1:
        #    pdb.set_trace()
        # re-key
        for group in temp_combined:
            if len(group) > 1:
                gs = group[0][2]
                ge = group[-1][2]
                gp = group[0][4]
                gf = group[-1][5]
            else:
                gs, ge = group[0][2], group[0][3]
                gp, gf = group[0][4], group[0][5]
            name = ''
            for pos,member in enumerate(group):
                if pos + 1 < len(group):
                    dist = group[pos + 1][3] - group[pos][3]
                    if dist > 1:
                        spacer = '...'
                    else:
                        spacer = ''
                else:
                    spacer = ''
                name += '%s(%s)%s' % (member[0], (member[3]-member[2])/len(member[0]), spacer)
            record.combined[name] = (((gs, ge), gp, gf),)
    return record

def worker(id, record, motifs, conf):
    # find msats
    record = msatSearch(record, motifs)
    # mask msat sequence
    record = softmask(record)
    # combined loci
    if conf.getboolean('GeneralParameters', 'CombineLoci'):
        record = combineLoci(record, conf.getint('GeneralParameters', 'CombineLociDist'))
        #pdb.set_trace()
    # open connection to microsatellites table
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    # add new data for mask table (which = microsatellite repeats)
    for match in record.matches:
        for repeat in record.matches[match]:
            motif_count = (repeat[0][1] - repeat[0][0])/len(match)
            #pdb.set_trace()
            cur.execute('''INSERT INTO mask (id, name, motif, start, end, 
            preceding, following, motif_count) VALUES (%s, %s,%s,%s,%s,%s,%s,%s) 
            ''', (id, record.id, match, repeat[0][0], repeat[0][1], repeat[1], 
            repeat[2], motif_count))
    if conf.getboolean('GeneralParameters', 'CombineLoci'):
        for combined in record.combined:
            for repeat in record.combined[combined]:
                cur.execute('''INSERT INTO combined (id, motif, start, end, 
                preceding, following) VALUES (%s, %s, %s, %s, %s, %s)
                ''', (id, combined, repeat[0][0], repeat[0][1], repeat[1], 
                repeat[2]))
    # update blob in 'main' table
    if record.matches:
        record_pickle = cPickle.dumps(record,1)
        # replace sequence with masked version in seq_trimmed - we can always
        # call UPPER() on it for non-masked version
        cur.execute('''UPDATE sequence SET seq_trimmed = %s, record = %s, msat = %s WHERE id = %s''', (record.seq.masked, record_pickle, True, id))
    else:
        cur.execute('''UPDATE sequence SET msat = %s WHERE id = %s''', (False, id))
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
    # TODO:  Switch index to reference sequence.id versus autoincrement value and create an index on it
    cur.execute('''CREATE TABLE mask (id INT UNSIGNED NOT NULL, 
        name VARCHAR(100), motif VARCHAR(8), start 
        MEDIUMINT UNSIGNED, end MEDIUMINT UNSIGNED, preceding MEDIUMINT 
        UNSIGNED, following MEDIUMINT UNSIGNED, motif_count SMALLINT 
        UNSIGNED, FOREIGN KEY (id) REFERENCES sequence (id), INDEX 
        mask_name (name)) ENGINE=InnoDB''')

def createCombinedLoci(cur):
    try:
        cur.execute('''DROP TABLE combined''')
    except:
        pass
    # TODO:  Switch index to reference sequence.id versus autoincrement value and create an index on it
    cur.execute('''CREATE TABLE combined (
        id INT UNSIGNED NOT NULL, 
        motif TEXT, 
        start MEDIUMINT UNSIGNED, 
        end MEDIUMINT UNSIGNED,
        preceding MEDIUMINT UNSIGNED,
        following MEDIUMINT UNSIGNED,
        FOREIGN KEY (id) REFERENCES sequence (id) 
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

def main():
    start_time = time.time()
    options, arg = interface()
    motd()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    # get sequence from dbase.  Use the server-side cursor to keep the load
    # down from large result sets
    #conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", 
    #db="454_msatcommander", cursorclass=MySQLdb.cursors.SSCursor)
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    createMaskTable(cur)
    updateSequenceTable(cur)
    if conf.getboolean('GeneralParameters', 'CombineLoci'):
        createCombinedLoci(cur)
    conn.commit()
    cur.execute('''SELECT id, record FROM sequence WHERE n_count <= 2 AND
    trimmed_len > 40''')
    #pdb.set_trace()
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
        print 'Multiprocessing.  Number of processors = %s\n' % n_procs
        # to test with fewer sequences
        #count = 0
        threads = []
        # access the data on sequence by sequence basis to avoid reading the 
        # entire table contents into memory
        rowcount = cur.rowcount
        pb = progress.bar(0,rowcount,60)
        pb_inc = 0
        row = cur.fetchone()
        while row is not None:
        #index = 0
        #while index < 5000:
            if len(threads) < n_procs:
                #row = cur.fetchone()
                # convert BLOB back to sequence record
                id, record = row[0], cPickle.loads(row[1])
                p = multiprocessing.Process(target=worker, args=(id, record, motifs, conf))
                p.start()
                threads.append(p)
                row = cur.fetchone()
                if (pb_inc+1)%1000 == 0:
                    pb.__call__(pb_inc+1)
                elif pb_inc + 1 == rowcount:
                    pb.__call__(pb_inc+1)
                pb_inc += 1
            else:
                for t in threads:
                    if not t.is_alive():
                        threads.remove(t)
    else:
        print 'Not using multiprocessing\n'
        # access the data on sequence by sequence basis to avoid 
        # reading the entire table contents into memory
        rowcount = cur.rowcount
        pb = progress.bar(0,rowcount,60)
        pb_inc = 0
        row = cur.fetchone()
        while row is not None:
            # convert BLOB back to sequence record
            id, record = row[0], cPickle.loads(row[1])
            #pdb.set_trace()
            worker(id, record, motifs, conf)
            row = cur.fetchone()
            if (pb_inc+1)%1000 == 0:
                pb.__call__(pb_inc+1)
            elif pb_inc + 1 == rowcount:
                pb.__call__(pb_inc+1)
            pb_inc += 1
    print '\n'
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'


if __name__ == '__main__':
    main()

