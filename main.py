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
import sqlite3
import ConfigParser
from Bio import Seq
from Bio import pairwise2
from Bio.SeqIO import QualityIO
from Bio.Alphabet import SingleLetterAlphabet


def revComp(seq):
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq.translate(bases)[::-1]

def revTags(tags):
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
        record = record[left.end():right.start()]
    elif left:
        record = record[left.end():]
    elif right:
        record = record[:right.start()]
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
        return high_score['tag'], high_score['matches'], high_score['seq_match'], high_score['seq_match_span'], high_score['start']
    else:
        return None, None, None, None, None

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
    # trim them off
    trimmed = trim(record, left_trim, right_trim)
    return trimmed

def midTrim(record, tags, max_gap_char=5, **kwargs):
    #if record.id == 'FX5ZTWB02DF2VB':
    #    pdb.set_trace()
    s = str(record.seq)
    mid, mid_trim = tagRegexer(s, tags)
    if mid_trim:
        trimmed = trim(record, mid_trim)
        return mid, trimmed, mid, 'Regex'
    elif kwargs['fuzzy']:
        mid, matches, seq_match, seq_match_span, start = smithWaterman(s, tags, 1)
        # we need to ensure that we aren't matching a potential mid-like
        # sequence in the latter portions of the individual sequence
        # so we check to make sure the MID start is within max_gap_char bases
        if matches and start < max_gap_char:
            #pdb.set_trace()
            tag_re = re.compile(('%s') % (seq_match_span))
            left_trim = re.search(tag_re, s)
            trimmed = trim(record, left_trim)
            return mid, trimmed, seq_match_span, 'SmithWaterman'
    return None, None, None, None

def linkerTrim(record, tags, max_gap_char=5, **kwargs):
    #if record.id == 'MID_NoError_SimpleX_OneError_FandRevcomp':
    #    pdb.set_trace()
    s = str(record.seq)
    link, left_link_trim, right_link_trim = tagRegexer(s, tags, False, False, True)
    # regex is only searching beginning and end of string.  i,e. 
    # "^Tag_Sequence" and "Tag_sequence$"
    if left_link_trim and right_link_trim:
        trimmed = trim(record, left_link_trim, right_link_trim)
        return link, trimmed, link, tags[link], 'RegexBoth'
    elif left_link_trim:
        trimmed = trim(record, left_link_trim)
        return link, trimmed, link, tags[link], 'RegexForward'
    elif right_link_trim:
        trimmed = trim(record, None, right_link_trim)
        return link, trimmed, link, tags[link], 'RegexReverseComp'
    elif kwargs['fuzzy']:
        l_tag, left_matches, left_seq_match, left_seq_match_span, left_start = smithWaterman(s, tags, 1)
        # check for reverse complement of linker
        reverseTags = revTags(tags)
        r_tag, right_matches, right_seq_match, right_seq_match_span, right_start = smithWaterman(s, reverseTags, 1)
        if (l_tag and r_tag) and (l_tag == r_tag):
            if left_start < max_gap_char and right_start > (len(s) - (len(r_tag) + max_gap_char)):
                left_link_re = re.compile(('%s') % (left_seq_match_span))
                right_link_re = re.compile(('%s') % (right_seq_match_span))
                left_trim = re.search(left_link_re, s)
                right_trim = re.search(right_link_re, s)
                trimmed = trim(record, left_trim, right_trim)
                return l_tag, trimmed, left_seq_match_span, tags[l_tag], 'SmithWatermanBoth'
        elif l_tag:
            #pdb.set_trace()
            if left_start < max_gap_char:
                left_link_re = re.compile(('%s') % (left_seq_match_span))
                left_trim = re.search(left_link_re, s)
                trimmed = trim(record, left_trim)
                return l_tag, trimmed, left_seq_match_span, tags[l_tag], 'SmithWatermanForward'
        elif r_tag:
            if right_start > (len(s) - (len(r_tag) + max_gap_char)):
                right_link_re = re.compile(('%s') % (right_seq_match_span))
                right_trim = re.search(right_link_re, s)
                trimmed = trim(record, None, right_trim)
                # make sure to return the forward sequence of the reverse_complement
                # so that our dictionary keys match up (keys are all forward)
                return revComp(r_tag), trimmed, right_seq_match_span, tags[revComp(r_tag)], 'SmithWatermanReverseComp'
    return None, None, None, None, None

def reverse(items):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)

def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    mid, reverse_mid = dict(conf.items('MID')), reverse(conf.items('MID'))
    linkers, reverse_linkers = dict(conf.items('Linker')), reverse(conf.items('Linker'))
    reverse_mid[None] = None
    reverse_linkers[None] = None
    clust = conf.items('Clusters')
    # build tag library 1X
    tags = tagLibrary(mid, linkers, clust)
    ###########
    # Database Stuff - to be moved
    ###########
    try:
        os.remove('screening_db.sqlite')
    except:
        pass
    conn = sqlite3.connect('screening_db.sqlite')
    cur = conn.cursor()
    cur.execute('''CREATE TABLE sequence (
    id integer primary key autoincrement,
    name text,
    mid text,
    mid_seq text,
    mid_match text,
    mid_method text,
    linker text,
    linker_seq text,
    linker_match text,
    linker_method text,
    cluster text,
    n_count integer,
    untrimmed_len integer)''')
    conn.commit()
    # for each sequence
    for record in QualityIO.PairedFastaQualIterator(
    open(conf.get('Input','sequence'), "rU"), 
    open(conf.get('Input','qual'), "rU")):
        #print record.id
        # convert low-scoring bases to 'N'
        untrimmed_len = len(record.seq)
        qual_trimmed = qualTrimming(record, 10)
        N_count = str(qual_trimmed.seq).count('N')
        #pdb.set_trace()
        # search on 5' (left) end for MID
        mid_seq, mid_trimmed, mid_match, mid_method = midTrim(qual_trimmed, tags, fuzzy=True)
        #pdb.set_trace()
        if mid_trimmed:
            # if MID, search for exact matches (for and revcomp) on Linker
            # provided no exact matches, use fuzzy matching (Smith-Waterman) +
            # error correction to find Linker
            linker_seq, linker_trimmed, linker_match, linker_cluster, linker_method = linkerTrim(mid_trimmed, tags[mid_seq], fuzzy=True)
        else:
            linker_seq, linker_trimmed, linker_match, linker_cluster, linker_method = None, None, None, None, None
        cur.execute("INSERT INTO SEQUENCE (name, mid, mid_seq, mid_match, mid_method, linker, linker_seq, linker_match, linker_method, cluster, n_count, untrimmed_len) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (record.id, reverse_mid[mid_seq], mid_seq, mid_match, mid_method, reverse_linkers[linker_seq], linker_seq, linker_match, linker_method, linker_cluster, N_count, untrimmed_len))
        #print record.id, reverse_mid[mid_seq], mid_match, mid_method, reverse_linkers[linker_seq], linker_match, linker_method, linker_cluster
        conn.commit()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'
    

if __name__ == '__main__':
    main()

