import sys
import numpy as np
import statsmodels.stats.multitest as mt
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
import rpy2
import gzip

mass = importr('MASS')
stats = importr('stats')
base = importr('base')
    
peakFile = sys.argv[1] # File of "bait" peaks to analyze (format: Chromosome <\t> peak-center-position)
tssFile =  sys.argv[2] # File of "prey" TSSs to analyze (format: Chromosome <\t> TSS-center-position)
conFile = sys.argv[3] # Contact file in juicer short format 

'''
in my case: 
peakfile = All dREG peaks file
tssFile = Promoters file (dREG peaks less than 100 bps from a GeneCode TSS)
conFile = Contact file in juicer short format
'''

out = sys.argv[4]


DIST = 300000 #300000 # Search window length
CAP = 2500 # Capture window around locus

f = robjects.r["exp"]

# Get contact data from file

loci, peaks, TSSs, contactProbabilities, pvalues = [],[],[],[],[]

print "Getting loci for contact analysis"
#
#'''
#The following loop fills up the list of loci surrounding peaks and adds empty list to their corresponding 
#positions in thelist of contacts and promoters to list those contacts and TSSs that falls within these loci
#'''
for locus in open(peakFile): # locus line format: <chr>\t<position = center of peak>
    locus = locus.strip().split()
    start = max(0, int(locus[1]) - DIST)
    stop = int(locus[1]) + DIST
    loci.append((locus[0],start,stop,int(locus[1])))
    TSSs.append([])
    
# Get TSSs from file

print out[-6:-4], "Sorting Promoters"

#'''
#The following loop assign promoters to loci surrounding peaks so that each promoter
#is writen as its distance from the center of the peak, the chromosome and its position
#'''
ct = 0
for tss in open(tssFile): # TSS line format: <chr>\t<position = center of TSS peak>
    split = tss.strip().split()
    chrom, pos = split[0], int(split[1])
    for i in range(len(loci)):
        if chrom == loci[i][0]:
            if loci[i][1] < pos < loci[i][2]:
                distance = pos - loci[i][3]
                TSSs[i].append((distance,chrom,pos))
    ct += 1
    if i%100 == 0: 
        print ct

print out[-6:-4], "Scanning",len(loci),"loci for interactions!"





# Calculate observed and expected distributions for each locus 
print 'Analyzing', len(loci), 'loci'
counter = 0
for locus in open(peakFile): # locus line format: <chr>\t<position = center of peak>
    expect, plus, minus = [],[],[]
    locus = locus.strip().split()
    chrom, position = locus[0], int(locus[1])
    chrom = chrom.strip("chr")
    
    # read the contacts around loci from the splited temp_contact_file
    try:
        with gzip.open(conFile+"_IN_"+peakFile+"_folder/temp_"+str(counter)+".gz") as contact_fp:
        #with open(conFile+"_IN_"+peakFile+"_folder/temp_"+str(counter)) as contact_fp:
            contacts= np.array([l.strip().split() for l in contact_fp.readlines()])
            print "contacts with locus number", contacts[0][0]
            contacts = contacts[:,1:4]
            # Calculate parameters for expected distribution
            distance = contacts[:,2].astype(int) - contacts[:,1].astype(int)
            expect += list(distance)
            #print len(expect)
    except IOError:
        expect.append(1)
        
    empDist = ivect(expect)
    robjects.r.assign("empDist",empDist)
    
    P = robjects.r("ecdf(empDist)") # Empirical Cumulative Distribution Function for all contacts' distances in the locus
    robjects.r.assign("P",P)
    
    # Get empirical distribution for locus (for the selected bait, rather than all contacts, within the locus)
    
    baitStart, baitStop = (position-CAP), (position+CAP)
    # Get all plus and minus diractions distances of interactions with positions within a 5kb (CAP*2) window from (+ and -) the peak's center
    for contact in contacts:
        if baitStart <= int(contact[1]) <= baitStop:
            distance = int(contact[2])-int(contact[1])
            plus.append(distance)
        if baitStart <= int(contact[2]) <= baitStop:
            distance = int(contact[2])-int(contact[1])
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
        print out[-6:-4], 'main loop', counter
    
out = open(sys.argv[4], "w")

corrected = mt.multipletests(pvalues, method='fdr_bh')
print 'writing'
p_count = 0
for prob in contactProbabilities:
    outStr = str('chr' + str(prob[0])) +"\t"+ 'distal' + "\t" + str(prob[1]) +"\t"+ 'proximal' +"\t"+ str(prob[3]) +"\t"+ str(prob[4]) +"\t"+ str(prob[5]) +"\t"+ str(prob[6]) +"\t"+ str(prob[7]) +"\t"+ str(corrected[1][p_count]) +"\n"
    out.write(outStr)
    p_count += 1

