#!/usr/bin/env python
# encoding: utf-8
"""
cleaner.py

Created by Brant Faircloth on 2009-12-17.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import os, sys, pdb, time, MySQLdb, subprocess, ConfigParser, tempfile, cPickle, multiprocessing, textwrap, optparse  


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


def main():
    start_time      = time.time()
    options, arg    = interface()
    motd()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
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
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    # get the cluster of critters in the dbase
    # determine some maximum match count for each cluster (i,e. 5; meaning
    # that no reads matches >= 5 other reads)
    # for each cluster determine the goods reads
        ## 0-matches
            ### put in good reads table
        ## 1-matches
            ### determine which is longest
            ### put longest in good reads
            ### put other in bad reads
        ## 2-matches
            ### determine which is longest
            ### put longest in good reads
            ### put other in bad reads
        ## n-matches
            ### determine which is longest
            ### put longest in good reads
            ### put other in bad reads            
    
    

if __name__ == '__main__':
    main()

