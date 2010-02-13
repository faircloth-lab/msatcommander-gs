#!/usr/bin/env python
# encoding: utf-8
"""
assembly.py

Created by Brant Faircloth on 2009-07-14.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os
import pdb
import shutil
import sqlite3
import tempfile
import subprocess
from modifications import Ace


class Cluster:
    def __init__(self, row):
        # set the data parameter to None
        self.data = None
        self.base_id, self.base_name, self.base_seq = row
        self.flagged = False
    
    def input():
        '''generic sequence input not from database'''
        pass

    def query(self, conn, cur):
        '''collect the sequence cluster data'''
        self.conn = conn
        # create new cursor instance to we don't overwrite the old cursor instance
        self.cur = cur
        # check sequence flagging - this query needs to be run "live" at every 
        # iteration to get updated values from dbase (i,e. since we are actively 
        # updating the records below, using the values from the input row won't do)
        flagged = self.cur.execute("SELECT sequence.assemble_flag FROM \
        sequence WHERE sequence.id = ?", (self.base_id,)).fetchall()
        if flagged[0][0]:
            self.flagged = True
        else:
            self.data = self.cur.execute("select blat.t_name, sequence.seq_seq, \
            sequence.assemble_flag, sequence.assemble_skip, t_start, t_end from \
            blat, sequence where blat.t_name = sequence.seq_id and \
            blat.id = (?) order by score DESC", (self.base_id,)).fetchall()
            # return number of sequences in cluster
            self.len = len(self.data)
            # if we don't get any results, we have no data
            if self.len < 1:self.data = None
        
    def filter(self):
        '''filter out those sequence with sequence assemble_flag = True and
        increment their assemble_skip'''
        assert self.data, "There are no data for filtering"
        self.cluster = [(self.base_name,self.base_seq)]
        self.excluded = ()
        self.included = ((self.base_name,),)
        for q in self.data:
            #pdb.set_trace()
            name, seq, flag, skip = q[0:4]
            if flag == 'True':
                # dont add to cluster
                # increment skip counter
                skip += 1
                # record skipped seqs in object
                self.excluded += ((skip,name,),)
                # record other data???
            else:
                # add to cluster
                self.cluster.append((name,seq))
        # increment counts of any sequences that were skipped
        self.cur.executemany("UPDATE sequence SET assemble_skip = ? where \
        seq_id = ?", self.excluded)
        # we do not yet want to flag either the base sequence or any associated
        # sequences as being assembled into a contig, because we have no idea
        # (yet) which sequences will be in the final contig with the base
    
    def cap3(self):
        '''assemble'''
        tdir = tempfile.mkdtemp(prefix='msatcommander-454-')
        tf = tempfile.mkstemp(prefix='align-%s-' % self.base_name, suffix='.tmp', dir=tdir)
        tf_handle = open(tf[1], 'w')
        for seq_pair in self.cluster:
            tf_handle.write(('>%s\n%s\n') % (seq_pair[0],seq_pair[1]))
        tf_handle.close()
        cap3_result, cap3_error = subprocess.Popen(['cap3 %s' % tf[1]], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if not cap3_error:
            # move the ace file to an appropriate location
            # using shutil in case filesystems are diff.
            self.temp_acefile = '%s.cap.ace' % tf[1]
            self.cap3 = True
        else:
            self.cap3 = False
        # delete the temp directory - use rmtree so we don't have to iterate
        # over damn files and os.remove()
        #shutil.rmtree(tdir)
        #
        #self.acefile = os.path.join(acefile, destination)
        # mark the base sequence as succesfully assembling

    def baseContig(self, destination):
        '''
        
        find the primary contig containing the base seq, ignore the others
        
        '''
        self.included = ()
        # parse the contig files to reduce memory consumption
        contigs = Ace.parse(open(self.temp_acefile))
        # just set base_contig to None, so we can catch any assemblies that 
        # don't include it
        base_contig = None
        for c in contigs:
            for r in c.reads:
                if r.rd.name == self.base_name:
                    base_contig = c
                    break
        if not base_contig:
            self.aligned = False
        else:
            self.aligned = True
            # reset the contig name to the base_name
            base_contig.name = self.base_name
            # get the sequence ids in the base contig
            for r in base_contig.reads:
                self.included += ((self.base_name, r.rd.name,),)
            # turn the contig back into an ACEFileRecord
            ace_file = Ace.ACEFileRecord()
            ace_file.contigs.append(base_contig)
            ace_file.ncontigs = len(ace_file.contigs)
            ace_file.nreads = base_contig.nreads
            # write the ACEFile to destination
            Ace.write(ace_file, destination)
            # delete the temporary directory
            shutil.rmtree(os.path.dirname(self.temp_acefile))
    
    def flag(self):
        '''flag out sequences that have already been assembled into contigs 
        (sequences cannot be in > 1 contig)'''
        # mark sequence as having been assembled.  Sequence can only be in 1 contig
        self.cur.executemany("UPDATE sequence SET assemble_flag = 'True', \
        contig_name = (?) where seq_id = (?)", self.included)
        self.conn.commit()
    
    def fewMatches(self):
        '''figure out what to do with contigs with <=1 match (filter out 0's? which 
        should probably be screened out in the quality/linker trimming)'''
        pass

if __name__ == '__main__':
    pass

