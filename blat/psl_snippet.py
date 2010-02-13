def pslPercentId(psl, protein=False, mRNA=True):
    '''convert the percent sequence ID of an alignment from a single line of a 
    parsed PSL file.  Code adapted from 
    http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    millibad = 0
    if protein:
        sizeMul = 3
    else:
        sizeMul = 1
    qAliSize = sizeMul * (psl['qEnd'] - psl['qStart'])
    tAliSize = psl['tEnd'] - psl['tStart']
    aliSize = min(qAliSize, tAliSize)
    if aliSize <= 0:return 0
    else:
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
            if mRNA:
                sizeDif = 0;
            else:
                sizeDif = -sizeDif
        insertFactor = psl['qNumInsert']
        if not mRNA:
            insertFactor += psl['tNumInsert']
        total = (sizeMul * (psl['match'] + psl['repMatch'] + psl['misMatch']))
        if total != 0:
            milliBad = (1000 * (psl['misMatch']*sizeMul + insertFactor + round(3*log(1+sizeDif)))) / total
        percent = round(100 - milliBad * 0.1,0)
        return percent

def pslScore(psl, sizeMul = 1):
    '''convert the score of an alignment from a single line of a 
    parsed PSL file.  Code adapted from 
    http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    return sizeMul * (psl['match'] + (psl['repMatch'] >> 1)) - sizeMul * psl['misMatch'] - psl['qNumInsert'] - psl['tNumInsert']
    
names = ['match', 'misMatch', 'repMatch', 'Ns', 'qNumInsert', 'qGapBases', 'tNumInsert', 'tGapBases', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd']
#values = [265,1,0,0,1,8,2,3,'+','FX5ZTWB02D840A',274,0,274,'chr1',197195432,21937230,21937499]
#psl = dict(zip(names, values))

values = [56,1,0,0,0,0,1,2,'-','FX5ZTWB02D840A',274,126,183,'chr14',125194864,39916111,39916170]
values =[46,1,0,0,1,2,2,6,'+','FX5ZTWB02D840A',274,145,194,'chr4',155630120,140050042,140050095,4]
values =[19,0,0,0,1,1,0,0,'-','FX5ZTWB02DA7SM',282,148,168,'chr9',124076172,70845115,70845134]
values =[167,1,0,0,0,0,0,0,'+','FX5ZTWB02C9U01',197,0,168,'chr2',181748087,161610225,161610393]
values =[160,1,0,0,2,3,2,2,'-','FX5ZTWB02EGX4Q',167,3,167,'chr3',159599783,95628043,95628206]
values =[213,0,0,0,3,3,1,1,'+','FX5ZTWB02DE0S2',218,0,216,'chr4',155630120,77471302,77471516]
values=[11,0,0,0,0,0,0,0,'+','FX5ZTWB02DE0S2',218,64,75,'chr6',149517037,98184812,98184823]
values=[27,0,0,0,1,2,0,0,'+','FX5ZTWB02DV4QF',127,7,36,'chr18',90772031,78081185,78081212]
values=[62,0,0,0,1,7,0,0,'-','FX5ZTWB02DV4QF',127,0,69,'chr9',124076172,8274373,8274435]
values=[197,0,0,0,0,0,0,0,'-','FX5ZTWB02DDQI0',205,0,197,'chr3',159599783,150775968,150776165]
values=[118,0,0,0,2,2,0,0,'-','FX5ZTWB02EDINL',133,0,120,'chr8',131738871,118559496,118559614]
values=[55,0,0,0,1,1,0,0,'-','FX5ZTWB02DCWY8',149,0,56,'chr17',95272651,79964808,79964863]
values=[92,0,0,0,1,1,2,3,'-','FX5ZTWB02DCWY8',149,56,149,'chr19',61342430,31170743,31170838]
values=[34,0,0,0,0,0,0,0,'-','FX5ZTWB02EQDJ2',225,191,225,'chr1',197195432,160025582,160025616]
