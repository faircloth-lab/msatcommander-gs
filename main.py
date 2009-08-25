#!/usr/bin/env python
# encoding: utf-8
"""
parser.py

Created by Brant Faircloth on 2009-08-14.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import re
import os
import pdb
import time
import numpy
import string
#import sqlite3
import MySQLdb
import ConfigParser
import multiprocessing
from Bio import Seq
from Bio import pairwise2
from Bio.SeqIO import QualityIO
from Bio.Alphabet import SingleLetterAlphabet


def revComp(seq):
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq.translate(bases)[::-1]

def revCompTags(tags):
    revTags = {}
    for tag in tags:
        revTags[revComp(tag)] = tags[tag]
    return revTags

def tagLibrary(mids, linkers, clust):
    tl = {}
    for c in clust:
        m,l = c[0].replace(' ','').split(',')
        org = c[1]
        if mids[m] not in tl.keys():
            tl[mids[m]] = {linkers[l]:org}
        else:
            tl[mids[m]][linkers[l]] = org
    return tl
            
def trim(record, left=None, right=None):
    '''takes regular expression objects'''
    if left and right:
        record = record[left:right]
    elif left:
        record = record[left:]
    elif right:
        record = record[:right]
    return record

def matches(seq_match_span, tag_match_span, allowed_errors):
    # deal with case where tag match might be perfect, but extremely gappy, 
    #e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or seq_match_span.count('-') > allowed_errors:
        return 0, 0
    # we don't really want to penalize mismatches/gaps with -1, rather we
    # just don't want to count them as matches (=0).  The case above should
    # help with long extensions across gaps
    else:
        #pdb.set_trace()
        seq_array, tag_array = numpy.array(list(seq_match_span)), numpy.array(list(tag_match_span))
        matches = sum(seq_array == tag_array)
        error = sum(seq_array != tag_array)
        #return sum((1 if s == tag_match_span[i] else 0) for i, s in enumerate(seq_match_span))
        return matches, error

def smithWaterman(seq, tags, allowed_errors):
    '''Borrowed & heavily modified from http://github.com/chapmanb/bcbb/tree/master'''
    #if seq == 'AAAACGTACGTGCGGATCTCCCCCTCAGCCTTTCCCTTCTACCAGACTAAAGAGATAGATAGATAAGATAGTA':
    #    pdb.set_trace()
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, 
            tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
        seq_match_span, tag_match_span = seq_match[start:end], tag_match[start:end]
        match, errors = matches(seq_match_span, tag_match_span, allowed_errors)
        #errors = len(tag) - match
        if match >= len(tag)-allowed_errors and match > high_score['matches'] and errors <= high_score['errors']:
            high_score['tag'] = tag
            high_score['seq_match'] = seq_match
            high_score['tag_match'] = tag_match
            high_score['score'] = score
            high_score['start'] = start
            high_score['end'] = end
            high_score['matches'] = match
            high_score['seq_match_span'] = seq_match_span
            high_score['errors'] = errors
    if high_score['matches']:
        return high_score['tag'], high_score['matches'], high_score['seq_match'], high_score['seq_match_span'], high_score['start'], high_score['end']
    else:
        return None

def tagRegexer(seq, tags, left=True, right=False, both=False):
    for tag in tags:
        if left:
            tag_re = re.compile(('^%s') % (tag))
            tag_trim = re.search(tag_re, seq)
            if tag_trim:
                return tag, tag_trim
        elif right:
            tag_re = re.compile(('%s$') % (revComp(tag)))
            tag_trim = re.search(tag_re, seq)
            if tag_trim:
                return tag, tag_trim
        elif both:
            left_tag_re = re.compile(('^%s') % (tag))
            right_tag_re = re.compile(('%s$') % (revComp(tag)))
            left_tag_trim = re.search(left_tag_re, seq)
            right_tag_trim = re.search(right_tag_re, seq)
            if left_tag_trim or right_tag_trim:
                return tag, left_tag_trim, right_tag_trim
    if left or right:
        return None, None
    else:
        return None, None, None

def qualTrimming(record, min_score=10):
    s = str(record.seq)
    sl = list(s)
    for q in enumerate(record.letter_annotations["phred_quality"]):
        if q[1] < min_score:
            sl[q[0]] = 'N'
    s = ''.join(sl)
    # find runs of ambiguous bases at 5' and 3' ends
    left_re, right_re = re.compile('^N+'),re.compile('N+$')
    left_trim, right_trim = re.search(left_re, s), re.search(right_re, s)
    if left_trim:
        left_trim = left_trim.end()
    if right_trim:
        right_trim = right_trim.end()
    return trim(record, left_trim, right_trim)

def midTrim(record, tags, max_gap_char=5, **kwargs):
    #if record.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    s = str(record.seq)
    mid = leftLinker(s, tags, max_gap_char, fuzzy=kwargs['fuzzy'])
    if mid:
        trimmed = trim(record, mid[3])
        tag, m_type, seq_match = mid[0], mid[1], mid[4]
        return tag, trimmed, seq_match, m_type
    else:
        return None

def SWMatchPos(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop

def leftLinker(s, tags, max_gap_char, **kwargs):
    for tag in tags:
        r = re.compile(('^%s') % (tag))
        match = re.search(r, s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    #if s == 'ACCTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG':
    #    pdb.set_trace()
    if not match and kwargs['fuzzy']:
        match = smithWaterman(s, tags, 1)
        # we can trim w/o regex
        if match:
            m_type = 'smithwaterman'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return tag, m_type, start, stop, seq_match
    else:
        return None

def rightLinker(s, tags, max_gap_char, **kwargs):
    #if s == 'GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG':
    #    pdb.set_trace()
    revtags = revCompTags(tags)
    for tag in revtags:
        r = re.compile(('%s$') % (tag))
        match = re.search(r, s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    if not match and kwargs['fuzzy']:
        match = smithWaterman(s, revtags, 1)
        # we can trim w/o regex
        if match:
            m_type = 'smithwaterman'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return revComp(tag), m_type, start, stop, seq_match
    else:
        return None

def linkerTrim(record, tags, max_gap_char=5, **kwargs):
    '''Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences'''
    #if record.id == 'MID_NoError_SimpleX_NoError_FandRevcomp':
    #    pdb.set_trace()
    m_type  = False
    s       = str(record.seq)
    left    = leftLinker(s, tags, max_gap_char=5, fuzzy=kwargs['fuzzy'])
    right   = rightLinker(s, tags, max_gap_char=5, fuzzy=kwargs['fuzzy'])
    if left and right and left[0] == right[0]:
        # we can have lots of conditional matches here
        if left[2] < max_gap_char and right[2] > (len(s) - (len(right[0]) + max_gap_char)):
            trimmed = trim(record, left[3], right[2])
            # left and right are identical so pass back the left info...
            tag, m_type, seq_match = left[0], left[1]+'-both', left[4]
        else:
            pass
    elif left and right and left[0] != right[0]:
        # flag
        pass
    elif left:
        if left[2] < max_gap_char:
            trimmed = trim(record, left[3])
            tag, m_type, seq_match = left[0], left[1], left[4]
        else:
            # flag
            pass
    elif right:
        if right[2] > (len(s) - (len(right[0]) + max_gap_char)):
            trimmed = trim(record, None, right[2])
            tag, m_type, seq_match = right[0], right[1], right[4]
        else:
            # flag
            pass
    if m_type:
        try:
            return tag, trimmed, seq_match, tags[tag], m_type
        except:
            pdb.set_trace()
    else:
        return None

def reverse(items):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)

def createSeqTable(c):
    try:
        c.execute('''DROP TABLE SEQUENCE''')
    except:
        pass
    c.execute('''CREATE TABLE SEQUENCE (id INT UNSIGNED NOT NULL AUTO_INCREMENT,name VARCHAR(100),mid VARCHAR(30),mid_seq VARCHAR(30),mid_match VARCHAR(30),mid_method VARCHAR(50),linker VARCHAR(30),linker_seq VARCHAR(30),linker_match VARCHAR(30),linker_method VARCHAR(50),cluster VARCHAR(50),n_count SMALLINT UNSIGNED,untrimmed_len SMALLINT UNSIGNED,PRIMARY KEY (id))''')

def worker(record, tags, reverse_mid, reverse_linkers):
    # we need a separate connection for each mysql cursor or they are going
    # start getting into locking hell and things go poorly. This is the
    # easiest solution.
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", db="454_msatcommander")
    cur = conn.cursor()
    #print os.getpid()
    # convert low-scoring bases to 'N'
    untrimmed_len = len(record.seq)
    qual_trimmed = qualTrimming(record, 10)
    N_count = str(qual_trimmed.seq).count('N')
    #pdb.set_trace()
    # search on 5' (left) end for MID
    mid = midTrim(qual_trimmed, tags, fuzzy=True)
    #pdb.set_trace()
    if mid:
        # if MID, search for exact matches (for and revcomp) on Linker
        # provided no exact matches, use fuzzy matching (Smith-Waterman) +
        # error correction to find Linker
        mid, trimmed, seq_match, m_type = mid
        #pdb.set_trace()
        linker = linkerTrim(trimmed, tags[mid], fuzzy=True)
        if linker:
            l_tag, l_trimmed, l_seq_match, l_critter, l_m_type = linker
        else:
            l_tag, l_trimmed, l_seq_match, l_critter, l_m_type = (None,) * 5
    else:
        mid, trimmed, seq_match, m_type = (None,) * 4
        l_tag, l_trimmed, l_seq_match, l_critter, l_m_type = (None,) * 5
    #pdb.set_trace()
    cur.execute("INSERT INTO SEQUENCE (name, mid, mid_seq, mid_match, mid_method, linker, linker_seq, linker_match, linker_method, cluster, n_count, untrimmed_len) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (record.id, reverse_mid[mid], mid, seq_match, m_type, reverse_linkers[l_tag], l_tag, l_seq_match, l_m_type, l_critter, N_count, untrimmed_len,))
    conn.commit()
    conn.close()
    return

def main():
    #pdb.set_trace()
    start_time = time.time()
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    mid, reverse_mid = dict(conf.items('MID')), reverse(conf.items('MID'))
    linkers, reverse_linkers = dict(conf.items('Linker')), reverse(conf.items('Linker'))
    reverse_mid[None] = None
    reverse_linkers[None] = None
    clust = conf.items('Clusters')
    # build tag library 1X
    tags = tagLibrary(mid, linkers, clust)
    # get multiprocessing
    n_procs = conf.get('Nprocs','processors')
    if n_procs == 'Auto':
        n_procs = multiprocessing.cpu_count() - 1
    else:
        n_procs = int(n_procs)
    print 'Multiprocessing.  Number of processors = ', n_procs
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    #create dbase table for records
    #connections = [MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", db="454_msatcommander") for i in range(n_procs)]
    # create a shitload of cursors - not sure if this will work
    #cursors = [i.cursor() for i in connections]
    #pdb.set_trace()
    # crank out a new table for the data
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", db="454_msatcommander")
    cur = conn.cursor()
    createSeqTable(cur)
    # for each sequence
    record = QualityIO.PairedFastaQualIterator(open(conf.get('Input','sequence'), "rU"), open(conf.get('Input','qual'), "rU"))
    mproc=False
    if mproc:
        count = 0
        try:
            jobs = [None] * n_procs
            while count < 5000:
                #pdb.set_trace()
                for i in range(n_procs):
                    count +=1
                    p = multiprocessing.Process(target=worker, args=(record.next(), tags, reverse_mid, reverse_linkers))
                    jobs[i]=p
                    p.start()
            for j in jobs:
                try:
                    j.join()
                except:
                    pass
        except StopIteration:
            pass
    else:
        count = 0
        try:
            while count < 5000:
                count +=1
                worker(record.next(), tags, reverse_mid, reverse_linkers)
        except StopIteration:
            pass
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'
    # let the threads finish writing to the dbase.
    #print "sleeping 10"
    #time.sleep(5)
    #for c in cursors:
    #    c.close()
    #db.close()

if __name__ == '__main__':
    main()

