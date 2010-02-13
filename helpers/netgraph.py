#!/usr/bin/env python
# encoding: utf-8
"""
netgraph.py

Created by Brant Faircloth on 2009-11-18.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import MySQLdb
import networkx as nx 
from networkx import graphviz_layout
import ConfigParser
from decimal import *


def nodetaker(nd, node):
    name, position = node.split(':')
    position = (int(position.split('-')[0]), int(position.split('-')[1]))
    if name not in nd.keys():
        nd[name] = [position]

def keychanger(d, q, k, zero, one):
    k_values = d[q][k]
    new_k = (zero, one)
    #new_k = [t_start, k[1]]
    del d[q][k]
    d[q][new_k] = k_values
    return d, new_k

def masternodehits(m):
    for k in m.keys():
        #pdb.set_trace()
        count = 0
        for i in m[k]:
            print k, i, len(m[k][i])

def masterslaves(m):
    for k in m.keys():
        #pdb.set_trace()
        count = 0
        for i in m[k]:
            print k,i
            for j in m[k][i]:
                print '\t', j

def overlaps((x1,x2), (y1,y2)):
    return not (x2 < y1 or y2 < x1)

def getseq(cur, seqid, (start, end)):
    cur.execute('''SELECT seq_trimmed from sequence_test where name = %s''', (seqid))
    sequence = cur.fetchone()
    return sequence[0][start:end]
    

def fastait(cur, seq, span, masternodes):
    print ">%s\n%s" % (seq, getseq(cur, seq,span))
    for s in masternodes[seq][span]:
        print '>%s\n%s' % (s[0], getseq(cur, s[0], (s[2],s[3])))


def nodeit(records, overlap = 10):
    #TODO:  check on double-ups within a particular masternode range???
    #TODO:  turn some of the shit below into functions?
    masternodes = {}
    slavenodes = {}
    for r in records:
        q_skip, t_skip, m_skip, entered = None, None, None, None
        q_name, t_name, percent, q_start, q_end, t_start, t_end = r[1],r[2],r[4],r[9],r[10],r[12],r[13]
        #pdb.set_trace()
        if percent >= Decimal("98.0"):
            #if q_name == 'FX5ZTWB04IPZ9T' and t_name in ['FX5ZTWB04I7TBS', 'FX5ZTWB04JM5Q4', 'FX5ZTWB04IVBC5', 'FX5ZTWB04JZE65']:  pdb.set_trace()
            # make sure the query node has not been used as a slave previously
            if q_name in slavenodes.keys():
                for k in slavenodes[q_name].keys():
                    if overlaps((q_start,q_end),k):
                        q_skip = True
                        break
            # make sure the target node has not been used as a slave previously
            if t_name in slavenodes.keys():
                for k in slavenodes[t_name].keys():
                    if overlaps((t_start, t_end), k):
                        t_skip = True
                        break
            # make sure the target node has not been used as a master previously
            if t_name in masternodes.keys():
                for k in masternodes[t_name].keys():
                    if overlaps((t_start, t_end), k):
                        m_skip = True
                        break
            if not q_skip and not t_skip and not m_skip:
                if q_name not in masternodes.keys():
                    # make sure there's no weird coverage/overlap issues
                    masternodes[q_name] = {(q_start,q_end):[(t_name, percent, t_start, t_end)]}
                else:
                    for k in masternodes[q_name].keys():
                        # check to see the overlap on the next record
                        # if it falls inside
                        if overlaps((q_start, q_end), k):
                            if k[1] >= q_end:
                                diff = q_end - k[0]
                            else:
                                diff = k[1] - q_start
                            if diff >= overlap:
                                masternodes[q_name][k].append((t_name, percent, t_start, t_end))
                                if q_start < k[0]:
                                    masternodes, k = keychanger(masternodes, q_name, k, q_start, k[1])
                                elif q_end > k[1]:
                                    masternodes, k = keychanger(masternodes, q_name, k, k[0], q_end)
                                entered = True
                                break     
                    if not entered:
                        masternodes[q_name][(q_start,q_end)]=[(t_name, percent, t_start, t_end)]
                    if t_name not in slavenodes.keys():
                        # stash the metadata of the slave node to compare against later
                        slavenodes[t_name] = {(t_start,t_end):(q_name, percent, q_start, q_end)}
                    else:
                        slavenodes[t_name][(t_start,t_end)]=(q_name, percent, q_start, q_end) 
    masternodehits(masternodes)
    #pdb.set_trace()
    return masternodes

def main():
    conf = ConfigParser.ConfigParser()
    conf.read('mc454.conf')
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    cur.execute('''SELECT * from graphtest2 order by percent DESC, length DESC''')
    records = cur.fetchall()
    masternodes = nodeit(records)
    pdb.set_trace()


if __name__ == '__main__':
    main()

