
#
#
# $Log$
#
#

import os
import sys
import random
import atexit

sys.path.append("/home/bbu-palin/nfsClustering/mabs-1_beta14/")

import align

# Parse commandline parameters
if len(sys.argv)<7:
    print """Usage: python ParameterOptimize.py Lambda xi mu nu nuc_per_rotation gff1 ..

    Give either positive initial values or negative fixed values"""
    sys.exit(0)

    

_system,hostname,vers,compTime,procType=os.uname()

# Number of alignments done for each bad and good
filesPerSide=6
GOROUNDS=60
limitNoImprove=6

# Initial settings
Lambda=2.0
xi=200.0
mu=0.2
nu=200.0
nuc_per_rotation=10.4


gLambda,gXi,gMu,gNu,gNuc=Lambda,xi,mu,nu,nuc_per_rotation


alignCount=0

#param=[Lambda,xi,mu,nu,nuc_per_rotation]
paramName=["Lambda","xi","mu","nu","nuc_per_rotation"]

paramMap={"Lambda":Lambda,"xi":xi,"mu":mu,"nu":nu,\
          "nuc_per_rotation":nuc_per_rotation}


OptimDims=[]
for i in range(1,6):
    try:
        val=float(sys.argv[i])
        paramMap[paramName[i-1]]=abs(val)
        if val<0.0:
            print >> sys.stderr,"Fixed %s=%g"%(paramName[i-1],abs(val))
        else:
            print >> sys.stderr,"Inital %s=%g"%(paramName[i-1],abs(val))
            OptimDims.append(paramName[i-1])
    except Exception,e:
#        print i,e
        raise
        

# Read the scanned matches from gffFiles

from cStringIO import StringIO

siteStrs={}
for fname in sys.argv[6:]:
    print "%s:"%(fname),
    siteStrs[fname]={}
    for seq,line in [(x.split()[0],x) for x in open(fname).readlines()]:
        if not siteStrs[fname].has_key(seq):
            print seq,
            siteStrs[fname][seq]=StringIO()
        siteStrs[fname][seq].write(line)
    print


# Select "good" files
gFiles=siteStrs.keys()
random.shuffle(gFiles)
gFiles=gFiles[:filesPerSide]
gFileSeq=[y[0].getvalue()+y[1].getvalue() for y in [siteStrs[x].values() for x in gFiles ]]
print "Selected %d good files: %s"%(len(gFileSeq),str(gFiles))


print "\nSELECTING BAD FILES:"
# select "bad" files
bFiles=[x for x in siteStrs.keys() if x not in gFiles]
random.shuffle(bFiles)
bFiles=gFiles+bFiles
#bFiles=bFiles[:filesPerSide]

bFileSeq=[]

while len(bFileSeq)<filesPerSide and len(bFiles)>0:
    file1=bFiles.pop()
    seqs1=siteStrs[file1].keys()
    random.shuffle(seqs1)
    file2=bFiles.pop()

    assert(len(seqs1)==2)
    if seqs1[0] in siteStrs[file2].keys() or seqs1[1] in siteStrs[file2].keys():
        print seqs1,"in",siteStrs[file2].keys()
        continue
    else:
        seq1=seqs1[0]
        seq2=random.choice(siteStrs[file2].keys())
        print "%s:%s <-> %s:%s"%(file1,seq1,file2,seq2)
        bFileSeq.append(siteStrs[file1][seq1].getvalue()+siteStrs[file2][seq2].getvalue())

print "Selected %d bad files"%(len(bFileSeq))
if len(bFileSeq)<filesPerSide:
    print "Hoped to get %d files"%(filesPerSide)

if len(bFileSeq)==0:
    print "Can't proceed! Please give more varied sequences"
    sys.exit(0)

sys.stdout.flush()
sys.stderr.flush()

def optimize(goodGFFstrs,badGFFstrs):
    """Function to compute the optimization values.
    Input: seqSet[fname][seq] the gff file like objects of good tf scan files
    goodFiles names of the good fnames used in the optimization
    badFileSeq (fname1,seq1,fname2,seq2) pairs used as background distribution.

    Returns (edge,goodScores,time)
    edge is the sum of ratios between good and bad alignments
    goodScores is number of good alignments better than the best bad alignment
    time is the time used in aligning"""
    
    timing=0.0
    goodScores=0
    edge=0.0
    badScores=0.0
    badMAX=0.0
    for badGFFstr in badGFFstrs:
        bgali=align.aligndata(badGFFstr,1,\
                              paramMap["Lambda"],paramMap["xi"],\
                              paramMap["mu"],paramMap["nu"],\
                              paramMap["nuc_per_rotation"])
        try:
            bad=bgali.nextBest()
        except AttributeError:
            print len(badGFFstr)
            print >> sys.stderr,"bad >%s<"%(bgali)
            raise
        timing+=bgali.secs_to_align
        badScores+=bad[0].score
        badMAX=max(bad[0].score,badMAX)
        print >> sys.stdout,"bad: %g %gsec"%(bad[0].score,bgali.secs_to_align)
        del(bgali)

    badAVG=badScores/len(badGFFstrs)
    
    for goodGFFstr in goodGFFstrs:
        print "Running good!"
        goodali=align.aligndata(goodGFFstr,100,
                              paramMap["Lambda"],paramMap["xi"],\
                              paramMap["mu"],paramMap["nu"],\
                              paramMap["nuc_per_rotation"])
        timing+=goodali.secs_to_align

        while 1:
            good=goodali.nextBest()
            goodScore=good[0].score
            print >> sys.stdout,"good: %g %gsec"%(good[0].score,goodali.secs_to_align)

            goodRatio=goodScore/badAVG

            if goodScore<badMAX: #goodRatio<1.0:
                break
            goodScores=goodScores+1
            print >> sys.stderr, "%d (%g) %g > %g"%(goodScores,goodRatio,goodScore,badMAX)
            edge=edge+goodRatio
        sys.stdout.flush()
        sys.stderr.flush()

    return (edge,goodScores,timing)


import math

def randomVector(n):
    return [random.random() for i in range(n) ]
##    v=[]
##    for i in range(n):
##        v.append(random.random())
##    return v

def randomUnitVector(n):
    norm=0.0
    v=[]
    for i in range(n):
        x=random.random()-0.5
        v.append(x)
        norm+=x**2
    norm=math.sqrt(norm)
    return map(lambda y:y/norm,v)


# Xi Mu Nu Nucl_per_rotation
randRateLimitsMap={"Lambda":100.0,"xi":1000.0,"mu":10.0,"nu":1000.0,"nuc_per_rotation":100.0}
ratesMap={"Lambda":10.0,"xi":10.0,"mu":1.0,"nu":10.0,"nuc_per_rotation":3.0}


absBestEdge=0.0
absBestParam=paramMap.copy()
def reportBest():
    print >> sys.stderr, "Best found: %g %s"%(absBestEdge,str(absBestParam))

atexit.register(reportBest)


# Greedy hill climbing search for best parameters.
def greedySearch():
    global absBestEdge
    global absBestParam
    bestEdge=-1.0
    noImprove=0
    totTime=0.0


    for iters in range(GOROUNDS):
        print >> sys.stderr,"params:",str(paramMap)
        #print >> sys.stderr,"best: %g %s"%(bestEdge,str(absBestParam))
        #edge,count,timing=1,2,3
        edge,count,timing=optimize(gFileSeq,bFileSeq)
        totTime+=timing
        toGoTime=(totTime/(iters+1)*GOROUNDS-totTime)
        print >> sys.stderr, "This %dth edge calculation (time elapsed %dh%dm, to go %dh%dm): %g"%(iters,totTime/3600,(totTime%3600)/60,toGoTime/3600,(toGoTime%3600)/60,edge)

        if edge>bestEdge:  # If edge is improved, store new parameters
            bestParam=paramMap.copy()
            bestEdge=edge
            if bestEdge>absBestEdge:
                absBestParam=bestParam.copy()
                absBestEdge=bestEdge
            print >> sys.stderr,"new best: %g %s"%(bestEdge,str(paramMap))
            #rates=map(lambda x:max(math.sqrt(x),1.0),bestParam[bginParam:endParam])
            for d in OptimDims:
                ratesMap[d]=max(math.sqrt(bestParam[d]),1.0)
            noImprove=0
        else:
            noImprove+=1

        if noImprove>limitNoImprove:
            print >> sys.stderr,"No improvement for %d rounds. Starting over.\nFinal best parameters with edge %g: %s"%(noImprove,bestEdge,str(bestParam))


            # Take a random jump to somewhere in the search space
            for d,x in zip(OptimDims,randomVector(len(OptimDims))):
                bestParam[d]=randRateLimitsMap[d]*x
                ratesMap[d]=max(math.sqrt(bestParam[d]),1.0)
                
            bestEdge=0.0
            noImprove=0

        # Take a random step somewhere close
        for d,r in zip(OptimDims,randomUnitVector(len(OptimDims))):
            ChgDir=ratesMap[d]*r
            paramMap[d]=max(ChgDir+bestParam[d],0.0)

        sys.stdout.flush()
        sys.stderr.flush()

    return (absBestEdge,absBestParam)


optParam=greedySearch()
print >> sys.stderr,"Absolutely best parameters:",optParam



