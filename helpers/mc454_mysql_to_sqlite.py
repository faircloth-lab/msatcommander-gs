
import os
import oursql
import sqlite3
import optparse
import ConfigParser

import pdb

def interface():
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')

    p.add_option('--output-db', dest = 'out_db', action='store', \
type='string', default = None, help='The path to the output db.', \
metavar='FILE')
    
    (options,arg) = p.parse_args()
    options.conf = os.path.expanduser(options.conf)
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg


def create_sqlite_tables(db, tables):
    """docstring for create_sqlite_tables"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute('''PRAGMA foreign_keys = ON''')
    try:
        for table in tables:
            if table == "sequence":
                c.execute('''CREATE TABLE sequence (id integer, name text, primary key (id))''')
            elif table == "mask" or table == "msats":
                c.execute('''CREATE TABLE msats (sequence_id integer NOT NULL,id integer NOT
                    NULL,motif text,start integer,end integer,preceding integer,following integer,
                    motif_count integer,PRIMARY KEY (sequence_id,id), FOREIGN KEY (sequence_id)
                    REFERENCES sequence (id))''')
            elif table == "combined":
                c.execute('''CREATE TABLE combined (sequence_id integer NOT NULL,id integer NOT NULL, motif
                    text,start integer,end integer,preceding integer,following integer,members
                    integer,PRIMARY KEY (sequence_id,id),FOREIGN KEY (sequence_id) REFERENCES
                    sequence (id))''')
            elif table == "combined_components":
                c.execute('''CREATE TABLE combined_components (sequence_id integer NOT NULL,combined_id
                    integer NOT NULL,motif text,length integer NOT NULL,FOREIGN KEY (sequence_id,
                    combined_id) REFERENCES combined (sequence_id, id))''')
            elif table == "primers":
                c.execute('''CREATE TABLE primers (sequence_id integer NOT NULL,id integer NOT NULL,primer
                    integer NOT NULL,left_p text NOT NULL,left_sequence text NOT NULL,left_tm
                    float NULL,left_gc float NULL,left_self_end float NULL,left_self_any float
                    NULL,left_hairpin float NULL,left_end_stability float NULL,left_penalty float
                    NULL,right_p text NOT NULL,right_sequence text NOT NULL,right_tm float
                    NULL,right_gc float NULL,right_self_end float NULL,right_self_any float
                    NULL,right_hairpin float NULL,right_end_stability float NULL,right_penalty
                    float NULL,pair_product_size float NULL,pair_compl_end float
                    NULL,pair_compl_any float NULL,pair_penalty float NULL,FOREIGN KEY
                    (sequence_id, id) REFERENCES combined (sequence_id, id))''')
            elif table == "tagged_primers":
                c.execute('''CREATE TABLE tagged_primers (sequence_id integer NOT NULL,
                    id integer NOT NULL,primer integer NOT NULL,best integer NOT NULL,tag text
                    NULL,tagged text NULL,tag_seq text NULL,common text NULL,pigtail_tagged text
                    NULL,pigtail_tag_seq text NULL,pigtail_common text NULL,left_p text
                    NULL,left_sequence text NULL,left_self_end float NULL,left_self_any float
                    NULL,left_hairpin float NULL,left_penalty float NULL,right_p text
                    NULL,right_sequence text NULL,right_self_end float NULL,right_self_any float
                    NULL,right_hairpin float NULL,right_penalty float NULL,pair_product_size float
                    NULL,pair_compl_end float NULL,pair_compl_any float NULL,pair_penalty float
                    NULL,FOREIGN KEY (sequence_id, id) REFERENCES combined (sequence_id, id))''')
    except sqlite3.OperationalError, e:
        if 'already exists' in e[0]:
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_sqlite_tables(db, tables)
            else:
                sys.exit(2)
    conn.commit()
    return conn, c

def get_mysql_tables(c):
    c.execute("SHOW TABLES")
    tables = c.fetchall()
    ordered_tables = []
    for item in ["sequence", "mask", "combined", "combined_components", "primers", "tagged_primers"]:
        if (item,) in tables:
            ordered_tables.append(item)
    return ordered_tables

def get_mysql_results(c, table):
    query = "SELECT * FROM {0}".format(table)
    c.execute(query)
    return c.fetchall()

def insert_to_sqlite(c, table, rows):
    qmark_len = len(rows[0])
    qmarks = '?,' * qmark_len
    qmarks = qmarks.rstrip(',')
    if table == "mask":
        table = "msats"
    query = "INSERT INTO {0} VALUES ({1})".format(table, qmarks)
    for row in rows:
        c.execute(query, row)
    return

def main():
    options, arg = interface()
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    mysql_con = oursql.connect(
        user = conf.get('Database','USER'), 
        passwd = conf.get('Database','PASSWORD'), 
        db = conf.get('Database','DATABASE')
        )
    mysql_cur = mysql_con.cursor()
    mysql_tables = get_mysql_tables(mysql_cur)
    sqlite_con, sqlite_cur = create_sqlite_tables(options.out_db, mysql_tables)
    for table in mysql_tables:
        rows = get_mysql_results(mysql_cur, table)
        insert_to_sqlite(sqlite_cur, table, rows)
    sqlite_con.commit()
    sqlite_con.close()
    mysql_con.close()
    
    
if __name__ == '__main__':
    main()
