import sys
import statsmodels.stats.multitest as mt
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
import rpy2

mass = importr('MASS')
stats = importr('stats')
base = importr('base')

peakFile = sys.argv[1] # File of "bait" peaks to analyze (format: Chromosome <\t> peak-center-position)
tssFile = sys.argv[2] # File of "prey" TSSs to analyze (format: Chromosome <\t> TSS-center-position)
conFile = sys.argv[3] # Contact file in juicer short format 

'''
in my case: 
peakfile = All dREG peaks file
tssFile = Promoters file (dREG peaks less than 100 bps from a GeneCode TSS)
conFile = Contact file in juicer short format
'''

out = open(sys.argv[4], "w")


DIST = 300000 # Search window length
CAP = 2500 # Capture window around locus

f = robjects.r["exp"]

# Get contact data from file

loci, contacts, peaks, TSSs, conList, contactProbabilities, pvalues = [],[],[],[],[],[],[]

print "Getting loci for contact analysis"

'''
The following loop fills up the list of loci surrounding peaks and adds empty list to their corresponding 
positions in thelist of contacts and promoters to list those contacts and TSSs that falls within these loci
'''
for locus in open(peakFile): # locus line format: <chr>\t<position = center of peak>

    locus = locus.strip()
    locus = locus.split()
    if int(locus[1]) >= DIST:
        start = int(locus[1]) - DIST
    else:
        start = 0
    stop = int(locus[1]) + DIST
    loci.append((locus[0],start,stop,int(locus[1])))
    contacts.append([])
    TSSs.append([])
    
# Get TSSs from file

print str(sys.argv[4])[-6:-4], "Sorting Promoters"

'''
The following loop assign promoters to loci surrounding peaks so that each promoter
is writen as its distance from the center of the peak, the chromosome and its position
'''
ct = 0
for tss in open(tssFile): # TSS line format: <chr>\t<position = center of TSS peak>
    split = tss.strip()
    split = tss.split()
    chrom, pos = split[0], int(split[1])

    for i in range(len(loci)):
        if chrom == loci[i][0]:
            if loci[i][1] < pos < loci[i][2]:
                distance = pos - loci[i][3]
                TSSs[i].append((distance,chrom,pos))
    ct += 1
    if i%100 == 0: 
        print ct

print str(sys.argv[4])[-6:-4], "Scanning",len(loci),"loci for interactions!"

'''
The following loop writes all contacts into a list as such: side1chr, side1pos, side2chr, side2pos
'''
# Make sure that the Juicer file format matches the line format as described here (short format):
# <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>
rowNumber = 0
for line in open(conFile):

    rowNumber += 1
        
    if rowNumber%100000 == 0: 
        print str(sys.argv[4])[-6:-4], rowNumber

    split = line.split()

    chromo1 = split[1]
    if chromo1[:3] == 'chr': chromo1 = chromo1[3:]
    chromo2 = split[5]
    if chromo2[:3] == 'chr': chromo2 = chromo2[3:]
    if chromo1 == chromo2:
        if int(split[2]) <= int(split[6]):
            conList.append([chromo1,split[2],chromo2,split[6]])
        else:
            conList.append([chromo2,split[6],chromo1,split[2]])

print str(sys.argv[4])[-6:-4], 'scanning', len(conList), 'contacts, to keep only contacts around loci'
# only keep contacts around loci
'''
The following loop assign to every locus around a peak, all the contacts that falls within it # this takes days to do, the assignement
'''
ce = 0
for entry in conList: 
    
    chr1,pos1,chr2,pos2 = 'chr' + str(entry[0]),int(entry[1]),'chr' + str(entry[2]),int(entry[3])
    for i in range(len(loci)):
        if chr1 == loci[i][0] and (int(loci[i][1]) <= pos1 <= int(loci[i][2])):
            contacts[i].append([chr1,pos1,pos2])
        elif chr1 == loci[i][0] and (int(loci[i][1]) <= pos2 <= int(loci[i][2])):
            contacts[i].append([chr1,pos1,pos2])
    ce += 1
    if ce%100000 == 0: 
        print str(sys.argv[4])[-6:-4], 'contacts to loci', ce, 'out of', len(conList)
# Calculate observed and expected distributions for each locus 

print 'Analyzing', len(loci), 'loci'
counter = 0
for locus in open(peakFile): # locus line format: <chr>\t<position = center of peak>
    expect, plus, minus = [],[],[]
    locus = locus.strip().split()
    chrom, position = locus[0], int(locus[1])
    if chrom[:3] == 'chr':
        chrom = chrom[3:]
    # Calculate parameters for expected distribution
    for contact in contacts[counter]:
        distance = int(contact[2]) - int(contact[1])
        expect.append(distance)


    if len(expect) < 1:
        expect.append(1)

    empDist = ivect(expect)
    robjects.r.assign("empDist",empDist)
    
    P = robjects.r("ecdf(empDist)") # Empirical Cumulative Distribution Function for all contacts' distances in the locus
    robjects.r.assign("P",P)
    
    # Get empirical distribution for locus (for the selected bait, rather than all contacts, within the locus)
    
    baitStart, baitStop = (position-CAP), (position+CAP)
    # Get all plus and minus diractions distances of interactions with positions within a 5kb (CAP*2) window from (+ and -) the peak's center
    for contact in contacts[counter]:
        if baitStart <= int(contact[1]) <= baitStop:
            distance = contact[2]-contact[1]
            plus.append(distance)
        if baitStart <= int(contact[2]) <= baitStop:
            distance = contact[2]-contact[1]
            minus.append(distance)
            
    #print 'plus length is:',len(plus), 'minus length is:', len(minus)
        # test plus strand contacts with promoters

    for tss in TSSs[counter]:
        if tss[0] > 0: #making sure to exclude self-interactions
            #print "Testing plus!"
            tssChrom, tssPosition, tssDist = tss[1],tss[2], tss[0]
            preyStart, preyStop = (tssDist - CAP), (tssDist + CAP)
            #print tssPosition, preyStart, preyStop
            
            
            if (tssDist - (CAP*1.5)) > 0:
                expStart, expStop = (tssDist - (CAP*1.5)), (tssDist + (CAP*1.5))
            else:
                expStart, expStop = 0, (tssDist + (CAP*1.5))
            
            robjects.r.assign("expStart",expStart)
            robjects.r.assign("expStop", expStop)
            z = robjects.r("seq(expStart,expStop,by=1)")
            robjects.r.assign("z",z)
            p = robjects.r("P(z)")
            exp = p[-1]-p[0]
        
            observed = 0
            for distance in plus:
                if preyStart <= distance <= preyStop:
                    observed += 1

            
            expected = exp*len(plus)
            #print "Expected:", len(plus),str(exp),str(expected)

                
            #print 'observed:', observed, 'expected:', expected
            v = robjects.FloatVector([observed, expected, (len(plus)-observed), (len(plus)-expected)])
            robjects.r.assign("v", v)
            test = robjects.r("matrix(v, nrow = 2)")
            robjects.r.assign("matrix", test)
            #print robjects.r("matrix")
            try:
                robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
                result = robjects.r("exact")
                p_val = result[0][0]
            except:
                p_val = 999.99
                #print "Exception!"
            #print p_val
            contactProbabilities.append([chrom, position, tssChrom, tssPosition, tssDist, p_val, observed, expected])
            pvalues.append(p_val)
    # test minus strand contacts
    
    for tss in TSSs[counter]:
        if tss[0] < 0:
            #print "Testing minus!"
            tssChrom, tssPosition, tssDist = tss[1],tss[2], tss[0]
            preyStart, preyStop = (tssDist - CAP), (tssDist + CAP)
            #print tssPosition, preyStart, preyStop
            
            if (abs(tssDist) - (CAP*1.5)) > 0:
                expStart, expStop = (abs(tssDist) - (CAP*1.5)), (abs(tssDist) + (CAP*1.5))
            else:
                expStart, expStop = 0, (abs(tssDist) + (CAP*1.5))
            
            robjects.r.assign("expStart",expStart)
            robjects.r.assign("expStop", expStop)
            z = robjects.r("seq(expStart,expStop,by=1)")
            robjects.r.assign("z",z)
            p = robjects.r("P(z)")
            exp = p[-1]-p[0]
            
            observed = 0
            for distance in minus:
                if preyStart <= (-1*distance) <= preyStop:
                    observed += 1

            expected = exp*len(minus)
            #print "Expected:",len(minus),str(exp),str(expected)

            #print 'observed:', observed, 'expected:', expected
            v = robjects.FloatVector([observed, expected, (len(minus)-observed), (len(minus)-expected)])
            robjects.r.assign("v", v)
            test = robjects.r("matrix(v, nrow = 2)")
            robjects.r.assign("matrix", test)
            #print robjects.r("matrix")
            try:
                robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
                result = robjects.r("exact")
                p_val = result[0][0]
            except:
                p_val = 999.99
                #print "Exception!"

            #print p_val
            contactProbabilities.append([chrom, position, tssChrom, tssPosition, tssDist, p_val, observed, expected])
            pvalues.append(p_val)
        
    counter += 1
    if counter%100000 == 0:
        print str(sys.argv[4])[-6:-4], 'main loop', counter
    
corrected = mt.multipletests(pvalues, method='fdr_bh')
print 'writing', str(sys.argv[4])[-6:-4]
p_count = 0
for prob in contactProbabilities:
    outStr = str('chr' + str(prob[0])) +"\t"+ 'distal' + "\t" + str(prob[1]) +"\t"+ 'proximal' +"\t"+ str(prob[3]) +"\t"+ str(prob[4]) +"\t"+ str(prob[5]) +"\t"+ str(prob[6]) +"\t"+ str(prob[7]) +"\t"+ str(corrected[1][p_count]) +"\n"
    out.write(outStr)
    p_count += 1
