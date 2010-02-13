import time
from numpy import *


def str_array(s, qual):
    sa = array(list(s))
    qa = array(qual)
    qai = qa > 20
    ra = char.array(['N']*len(sa))
    w = putmask(ra, qai, sa)
    print ra

def enu(s, qual):
    sl = list(s)
    for q in enumerate(qual):
        if q[1] < 20:
            sl[q[0]] = 'N'
    print sl
	


qual = [37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 35, 33, 33, 32, 32, 27, 19, 15, 15, 15, 15, 14, 18, 13, 21, 21, 28, 28, 28, 30, 29, 29, 29, 30, 32, 33, 28, 27, 23, 23, 24, 27, 30, 30, 31, 29, 19, 19, 19, 19, 19, 19, 29, 26, 26, 22, 13, 13, 15, 23, 19, 24, 19, 19, 19, 32, 32, 32, 33, 33, 33, 32, 32, 31, 33, 33, 32, 23, 23, 23, 24, 19, 19, 19, 27, 31, 29, 27, 27, 27, 28, 29, 30, 29, 24, 19, 17, 17, 24]

s = 'TGTACTACTCTGATTGATAGATACATAGATACATAGATAGATACATAGATGGTAGATCAATAGATACATAGATAGATGCCTAGATAGATGCCTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGTAAAACTAGGCAGCTCCTGAGGAATGACACCTAAGGTTGCTTTCTGATTCACAGAGATCCTCCTGCTACTGCCTCCACGTGGA'

start = time.time()
str_array(s,qual)
end = time.time()
print 'array time: ', end - start
start = time.time()
enu(s,qual)
end = time.time()
print 'enu time: ', end - start
