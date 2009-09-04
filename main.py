#!/usr/bin/env python
# encoding: utf-8
"""
main.py

Created by Brant Faircloth on 2009-08-14.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import re, os, pdb, time, numpy, string, MySQLdb, ConfigParser, multiprocessing, cPickle
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

def allPossibleTags(mids, linkers, clust):
    at = []
    rat = []
    for c in clust:
        m,l = c[0].replace(' ','').split(',')
        at.append(linkers[l])
        rat.append(re.compile('%s' % linkers[l]))
        at.append(revComp(linkers[l]))
        rat.append(re.compile('%s' % revComp(linkers[l])))
    return at, rat
            
def trim(record, left=None, right=None):
    if left and right:
        record = record[left:right]
    elif left:
        record = record[left:]
    elif right:
        record = record[:right]
    return record

def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    # deal with case where tag match might be perfect, but extremely gappy, 
    #e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or seq_match_span.count('-')\
     > allowed_errors:
        return 0, 0
    else:
        #pdb.set_trace()
        seq_array, tag_array = numpy.array(list(seq_match_span)), \
        numpy.array(list(tag_match_span))
        matches = sum(seq_array == tag_array)
        error = sum(seq_array != tag_array) + (len(tag) - \
        len(tag_match_span.replace('-','')))
        # Original scoring method from 
        # http://github.com/chapmanb/bcbb/tree/master treats gaps incorrectly:
        #return sum((1 if s == tag_match_span[i] else 0) for i, s in 
        #enumerate(seq_match_span))
        return matches, error

def smithWaterman(seq, tags, allowed_errors):
    '''Borrowed & heavily modified from http://github.com/chapmanb/bcbb/tree/master'''
    #if seq == 'CGAGAGATACAAAAGCAGCAGCGGAATCGATTCCGCTGCTGC':
    #    pdb.set_trace()
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 
        'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, 
        tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
        seq_match_span, tag_match_span = seq_match[start:end], tag_match[start:end]
        match, errors = matches(tag, seq_match_span, tag_match_span, allowed_errors)
        if match >= len(tag)-allowed_errors and match > high_score['matches'] \
        and errors <= high_score['errors']:
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
        return high_score['tag'], high_score['matches'], \
        high_score['seq_match'], high_score['seq_match_span'], \
        high_score['start'], high_score['end']
    else:
        return None

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
    mid = leftLinker(s, tags, max_gap_char, True, fuzzy=kwargs['fuzzy'])
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

def leftLinker(s, tags, max_gap_char, gaps=False, **kwargs):
    for tag in tags:
        if gaps:
            r = re.compile(('^%s') % (tag))
        else:
            r = re.compile(('^[acgtnACGTN]{0,%s}%s') % (max_gap_char, tag))
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
            m_type = 'fuzzy'
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
        #r = re.compile(('%s$') % (tag))
        r = re.compile(('%s[acgtnACGTN]{0,%s}$') % (tag, max_gap_char))
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
            m_type = 'fuzzy'
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
    #if record.id == 'FX5ZTWB02DOPOT':
    #    pdb.set_trace()
    m_type  = False
    s       = str(record.seq)
    left    = leftLinker(s, tags, max_gap_char=5, fuzzy=kwargs['fuzzy'])
    right   = rightLinker(s, tags, max_gap_char=5, fuzzy=kwargs['fuzzy'])
    if left and right and left[0] == right[0]:
        # we can have lots of conditional matches here
        if left[2] < max_gap_char and right[2] > (len(s) - (len(right[0]) + \
        max_gap_char)):
            trimmed = trim(record, left[3], right[2])
            # left and right are identical so largely pass back the left
            # info... except for m_type which can be a combination
            tag, m_type, seq_match = left[0], left[1]+'-'+right[1]+'-both', \
            left[4]
        else:
            pass
    elif left and right and left[0] != right[0]:
        # flag
        if left[2] < max_gap_char and right[2] > (len(s) - (len(right[0]) + \
        max_gap_char)):
            trimmed = None
            tag, m_type, seq_match = None, 'tag-mismatch', None
    elif left:
        if left[2] < max_gap_char:
            trimmed = trim(record, left[3])
            tag, m_type, seq_match = left[0], left[1]+'-left', left[4]
        else:
            # flag
            pass
    elif right:
        if right[2] > (len(s) - (len(right[0]) + max_gap_char)):
            trimmed = trim(record, None, right[2])
            tag, m_type, seq_match = right[0], right[1]+'-right', right[4]
        else:
            # flag
            pass
    if m_type:
        try:
            return tag, trimmed, seq_match, tags[tag], m_type
        except:
            return tag, trimmed, seq_match, None, m_type
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
        c.execute('''DROP TABLE sequence_test''')
    except:
        pass
    c.execute('''CREATE TABLE sequence_test (id INT UNSIGNED NOT NULL 
        AUTO_INCREMENT,name VARCHAR(100),mid VARCHAR(30),mid_seq VARCHAR(30),
        mid_match VARCHAR(30),mid_method VARCHAR(50),linker VARCHAR(30),
        linker_seq VARCHAR(30),linker_match VARCHAR(30),linker_method 
        VARCHAR(50),cluster VARCHAR(50),concat_seq VARCHAR(30), 
        concat_match varchar(30), concat_method VARCHAR(50),
        n_count SMALLINT UNSIGNED, untrimmed_len SMALLINT UNSIGNED, 
        seq_trimmed TEXT, trimmed_len SMALLINT UNSIGNED, record BLOB, PRIMARY
        KEY (id))''')

def concatCheck(record, all_tags, all_tags_regex, reverse_linkers, **kwargs):
    s = str(record.seq)
    m_type = None
    # do either/or to try and keep speed up, somewhat
    #if not kwargs['fuzzy']:
    #pdb.set_trace()
    for tag in all_tags_regex:
        match = re.search(tag, s)
        if match:
            tag = tag.pattern
            m_type = 'regex-concat'
            seq_match = tag
            break
    if not match and ['fuzzy']:
    #else:
        match = smithWaterman(s, all_tags, 1)
        # we can trim w/o regex
        if match:
            tag = match[0]
            m_type = 'fuzzy-concat'
            seq_match = match[3]
    if m_type:
        return tag, m_type, seq_match
    else:
        return None, None, None
            
def worker(record, qual, tags, all_tags, all_tags_regex, reverse_mid, reverse_linkers):
    # we need a separate connection for each mysql cursor or they are going
    # start going into locking hell and things will go poorly. Creating a new 
    # connection for each worker process is the easiest/laziest solution.
    # Connection pooling (DB-API) didn't work so hot, but probably because 
    #I'm slightly retarded.
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", db="454_msatcommander")
    cur = conn.cursor()
    # convert low-scoring bases to 'N'
    untrimmed_len = len(record.seq)
    qual_trimmed = qualTrimming(record, qual)
    N_count = str(qual_trimmed.seq).count('N')
    # search on 5' (left) end for MID
    mid = midTrim(qual_trimmed, tags, fuzzy=True)
    #TODO:  Add length parameters
    if mid:
        # if MID, search for exact matches (for and revcomp) on Linker
        # provided no exact matches, use fuzzy matching (Smith-Waterman) +
        # error correction to find Linker
        mid, trimmed, seq_match, m_type = mid
        #TODO:  Add length parameters
        linker = linkerTrim(trimmed, tags[mid], fuzzy=True)
        if linker:
            l_tag, l_trimmed, l_seq_match, l_critter, l_m_type = linker
        else:
            l_tag, l_trimmed, l_seq_match, l_critter, l_m_type, concat_type, \
            concat_count = (None,) * 7
    else:
        mid, trimmed, seq_match, m_type = (None,) * 4
        l_tag, l_trimmed, l_seq_match, l_critter, l_m_type, concat_type, \
        concat_count = (None,) * 7
    # check for concatemers
    concat_check = False
    if concat_check:
        if l_trimmed and len(l_trimmed.seq) > 0:
            concat_tag, concat_type, concat_seq_match = concatCheck(l_trimmed, 
                all_tags, all_tags_regex, reverse_linkers, fuzzy=True)
        else:
            concat_tag, concat_type, concat_seq_match = None, None, None
    else:
        concat_tag, concat_type, concat_seq_match = None, None, None
    # if we are able to trim the linker
    if l_trimmed:
        record = l_trimmed
    # if we are able to trim the MID
    elif trimmed:
        record = trimmed
    # pickle the sequence record, so we can store it as a BLOB in MySQL, we
    # can thus recurrect it as a sequence object when we need it next.
    record_pickle = cPickle.dumps(record,1)
    cur.execute('''INSERT INTO sequence_test (name, mid, mid_seq, mid_match, 
        mid_method, linker, linker_seq, linker_match, linker_method, cluster, 
        concat_seq, concat_match, concat_method, n_count, untrimmed_len, 
        seq_trimmed, trimmed_len, record) 
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
        (record.id, reverse_mid[mid], mid, seq_match, m_type, 
        reverse_linkers[l_tag], l_tag, l_seq_match, l_m_type, l_critter, 
        concat_tag, concat_seq_match, concat_type, N_count, untrimmed_len,
        record.seq, len(record.seq), record_pickle))
    #pdb.set_trace()
    cur.close()
    conn.commit()
    # keep our connection load low
    conn.close()
    return

def main():
    start_time = time.time()
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    mid, reverse_mid = dict(conf.items('MID')), reverse(conf.items('MID'))
    linkers, reverse_linkers = dict(conf.items('Linker')), reverse(conf.items('Linker'))
    reverse_mid[None] = None
    reverse_linkers[None] = None
    clust = conf.items('Clusters')
    qual = conf.getint('Qual', 'MIN_SCORE')
    # build tag library 1X
    tags = tagLibrary(mid, linkers, clust)
    all_tags, all_tags_regex = allPossibleTags(mid, linkers, clust)
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    # crank out a new table for the data
    conn = MySQLdb.connect(user="python", passwd="BgDBYUTvmzA3", db="454_msatcommander")
    cur = conn.cursor()
    createSeqTable(cur)
    cur.close()
    conn.close()
    record = QualityIO.PairedFastaQualIterator(
    open(conf.get('Input','sequence'), "rU"), 
    open(conf.get('Input','qual'), "rU"))
    #pdb.set_trace()
    if conf.getboolean('Multiprocessing', 'MULTIPROCESSING'):
        # get num processors
        n_procs = conf.get('Multiprocessing','processors')
        if n_procs == 'Auto':
            # TODO:  change this?
            # we'll start 2X-1 threads (X = processors).
            n_procs = multiprocessing.cpu_count()
        else:
            n_procs = int(n_procs)
        print 'Multiprocessing.  Number of processors = ', n_procs
        # to test with fewer sequences
        #count = 0
        try:
            threads = []
            # this is approximately 150% faster than the code below.
            while record:
            #while count < 5000:
                if len(threads) < n_procs:
                    #count += 1
                    p = multiprocessing.Process(target=worker, args=(
                    record.next(), qual, tags, all_tags, all_tags_regex, 
                    reverse_mid, reverse_linkers))
                    p.start()
                    threads.append(p)
                else:
                    for t in threads:
                        if not t.is_alive():
                            threads.remove(t)
            ## TODO:  move this to more pool-based list of processes 
            ## (see example on intertubes)
            ##jobs = [None] * n_procs
            ## to test with fewer sequences
            ##while count < 5000:
            ## to test with all sequences
            ##while record:
            #    for i in range(n_procs):
            #        # to test with fewer sequences
            #        count +=1
            #        p = multiprocessing.Process(target=worker, args=(
            #        record.next(), qual, tags, all_tags, all_tags_regex, 
            #        reverse_mid, reverse_linkers))
            #        jobs[i]=p
            #        p.start()
            #    for j in jobs:
            #        try:
            #            j.join()
            #        except:
            #            pass
        except StopIteration:
            pass
    else:
        print 'Not using multiprocessing'
        #count = 0
        try:
            #while count < 1000:
            while record:
                #count +=1
                worker(record.next(), qual, tags, all_tags, all_tags_regex, 
                reverse_mid, reverse_linkers)
        except StopIteration:
            pass
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()

