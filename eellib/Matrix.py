# -*- coding: UTF-8 -*-
"""Python interface and logic to TFBS matrix matching"""
import operator
import string
import math

from eellib import _c_matrix

#
# $Log$
# Revision 1.15  2006/08/31 10:08:36  kpalin
# More Debug output and TF distance calculations.
#
# Revision 1.14  2005/05/19 07:49:35  kpalin
# Merged Waterman-Eggert style suboptimal alignments and
# SNP matching.
#
# Revision 1.13.2.1  2005/05/09 07:15:09  kpalin
# Matrix distances.
#
# Revision 1.13  2005/03/22 13:17:13  kpalin
# Merged some fixes surfacing from testing the public version.
#
# Revision 1.12  2005/03/08 11:19:49  kpalin
# Removed a debugging exception.
#
# Revision 1.11  2005/03/08 10:37:29  kpalin
# Fixed the use of markov Background, such that one background
# takes care of all Matrix-objects.
#
# Revision 1.10  2005/02/21 09:31:33  kpalin
# Tools to output matrices in alternative formats.
#
# Revision 1.9  2005/01/27 09:12:05  kpalin
# Added functions to format the matrix in string good
# for AHAB and as list of string scoring more than given limit
#
# Revision 1.8  2005/01/12 13:34:55  kpalin
# Added Tkinter/Tix Graphical user interface and command -no-gui to
# avoid it.
#
# Revision 1.7  2005/01/07 13:41:25  kpalin
# Works with py2exe. (windows executables)
#
# Revision 1.6  2004/12/22 11:14:24  kpalin
# Some fixes for better distributability
#
# Revision 1.5  2004/12/22 08:02:59  kpalin
# Hopefully more IO efficient TFBS search.
#
# Revision 1.4  2004/04/08 13:11:58  kpalin
# InfoContent and CVS loging.
#
#

from cStringIO import StringIO

def log2(x):
    "Base 2 logarithm"
    #return math.log(x)
    return math.log(x)/math.log(2.0)


class Matrix:
    "represents a binding site matrix"

    backGround=None
    freqA,freqC,freqG,freqT=0.25,0.25,0.25,0.25
    pseudoCount=1.0
    def __init__(self,filename):
        "reads matrix from file"

        File=open(filename,'r')
        self.fname=filename
        self.name=filename
        self.LLMatrix=filter(lambda x:len(x),[[string.atoi(entry)
                        for entry in line.split()]
                       for line in File.read().split('\n')])
        File.close()
        self.initWeights()
        #self.draw()


    def toMatCompare(self,s=None):
        "Return the matrix as file like object suitable for MatCompare input"
        if not s:
            s=StringIO()
        
        s.write("M.%s\n"%(self.name))
        self.toColumnMatrix(s)
        s.write("END")

        return s

    
    def toColumnMatrix(self,s=None):
        "Return a file like object having the transposed matrix"
        if not s:
            s=StringIO()
        
        for i in range(len(self)):
            A,C,G,T=self.LLMatrix[0][i],self.LLMatrix[1][i],self.LLMatrix[2][i],self.LLMatrix[3][i]
            s.write("%d\t%d\t%d\t%d\n"%(A,C,G,T))

        return s
        
    def toCisevolver(self,s=None):
        "Return the matrix as file like object formatted for cisevolver"
        if not s:
            s=StringIO()

        s.write("\n".join(["%f %f %f %f"%(x[0],x[1],x[2],x[3]) for x in zip(*self.freq)]))
        return s
    
    def toAhab(self,s=None):
        "Return the matrix as file like object formatted for ahab"
        if not s:
            s=StringIO()
        
        s.write(">%s Matrix\t%d\n"%(self.name,len(self)))
        self.toColumnMatrix(s)
        s.write("<")
        return s

    def seqsBetterThan(self,limit):
        "Return sequences scoring better than limit"
        if Matrix.backGround:
            raise ValueError("Can't use markov background")
        self.initWeights()
        ret=[] # Found sequences better than limit

        acgtWeights=[zip(x,['A','C','G','T']) for x in zip(*self.M_weight)]

        maxToGo=[]
        score=0.0
        acgtWeights.reverse()
        for i in acgtWeights:
            i.sort()
            i.reverse()
            maxToGo.append(score)
            score+=i[0][0]
        acgtWeights.reverse()
        maxToGo.reverse()

        if score<limit:
            return []
        
        def recurse(i=0,s=[],score=0.0):
            if i>=len(self):
                if score>=limit:
                    #print "".join(s)
                    ret.append(("".join(s),score))
                return
            acgt= acgtWeights[i]
            for sc,nucl in acgt:
                if score+sc+maxToGo[i]>limit:
                    recurse(i+1,s+[nucl],score+sc)
                else:
                    return
                
        recurse()
        return ret
        



    def getName(self):
        """Gives the name of this matrix"""
        return self.name

    def __len__(self):
        """Return the number of columns in this matrix"""
        return len(self.LLMatrix[0])

    def draw(self):
        "draws the matrix"
        for line in self.LLMatrix:
            for entry in line:
                print '%3.0f' % entry,
            print

    def drawWeights(self):
        "draws the matrix"
        for line in self.M_weight:
            for entry in line:
                print '%6.2f' % entry,
            print
            

    def setBG(self, seq):
        """Sets the 0-order background frequences from and example sequence."""
        countA,countC,countG,countT=seq.count("A")+seq.count("a"),seq.count("C")+seq.count("c"),seq.count("G")+seq.count("g"),seq.count("T")+seq.count("t")
        overall=countA+countC+countG+countT+0.0    # what happens with 'N'?
        self.setBGfreq(countA/overall,countC/overall,countG/overall,countT/overall)
        #self.freqA=countA/overall
        #self.freqC=countC/overall
        #self.freqG=countG/overall
        #self.freqT=countT/overall
        #self.freqX now contains the background self.frequency of base X
        self.initWeights()


    def initWeights(self):
        """Helper to initialize the matrix weights for 0- or higher order background models"""
        self.M_weight=[]
        self.computeInfoContent()
        self.trivialWeights()
        if Matrix.backGround:
            self.M_weight=[]
            self.positiveWeights()


    def positiveWeights(self):
        """Update weights for higher order markov background.

        Update weights only for positive probability.
        Background is taken care elsewhere"""

        sum=reduce(lambda s,row:map(operator.add,s,row),self.LLMatrix,[Matrix.pseudoCount]*len(self.LLMatrix[0]))


        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqA*Matrix.pseudoCount))/tot),
                            self.LLMatrix[0],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*Matrix.pseudoCount))/tot),
                            self.LLMatrix[1],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*Matrix.pseudoCount))/tot),
                            self.LLMatrix[2],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqT*Matrix.pseudoCount))/tot),
                            self.LLMatrix[3],sum))

        # Don't know how to compute maxscore with markov background. It varies.
        #maxscore is the highest reachable score (disregarding BG)
        #self.maxscore=0.0
        #for i in range(len(self.M_weight[0])):
        #    self.maxscore-=max(self.M_weight[0][i],
        #                         self.M_weight[1][i],
        #                         self.M_weight[2][i],
        #                         self.M_weight[3][i])

        



    def computeFreqMatrix(self):
        "Computes the self.freq matrix. The frequency table for this matrix."
        # Sum of the columns.
        sum=reduce(lambda s,row:map(operator.add,s,row),self.LLMatrix,[Matrix.pseudoCount]*len(self.LLMatrix[0]))


        # Frequencies of a nucleotide
        self.freq=[ [ (x+bgFreq*Matrix.pseudoCount)/(tot) for x,tot in zip(matrix,sum) ] \
                    for matrix,bgFreq in zip(self.LLMatrix,(self.freqA,self.freqC,self.freqG,self.freqT)) ]
        
    def computeInfoContent(self):
        "Computes and sets the information content for this matrix."

        self.computeFreqMatrix()

        # Compute information content of the motif
        InfoContent=[[ f*log2(f/bgF) for f in fLine] for (fLine,bgF) in zip(self.freq,(self.freqA,self.freqC,self.freqG,self.freqT)) ]


        self.InfoContent=reduce(operator.add,reduce(lambda x,y:x+y,InfoContent,[]))
    def trivialWeights(self):
        """Update weights according to 0-order background"""

        sum=reduce(lambda s,row:map(operator.add,s,row),self.LLMatrix,[Matrix.pseudoCount]*len(self.LLMatrix[0]))

        # Compute weights.
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqA*Matrix.pseudoCount))/
                                               (tot*self.freqA)),
                            self.LLMatrix[0],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*Matrix.pseudoCount))/
                                               (tot*self.freqC)),
                            self.LLMatrix[1],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqG*Matrix.pseudoCount))/
                                               (tot*self.freqG)),
                            self.LLMatrix[2],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqT*Matrix.pseudoCount))/
                                               (tot*self.freqT)),
                            self.LLMatrix[3],sum))


        #maxscore is the highest reachable score (which respect to the BG)
        self.maxscore=0.0
        for i in range(len(self.M_weight[0])):
            self.maxscore+=max(self.M_weight[0][i],
                                 self.M_weight[1][i],
                                 self.M_weight[2][i],
                                 self.M_weight[3][i])



    def minimumKLdistance(self,other):
        "Compute the minimum Kullback-Leibler distance between self and other"
        if len(self)<len(other):
            return other.minimumKLdistance(self)

        minKLdist=1e50
        self.computeFreqMatrix()
        other.computeFreqMatrix()

        assert(len(self)>=len(other))

        for shift in range(1-len(other),len(self)):
            kld=self.KLdistance(other,shift)
            #print shift,kld
            if kld<minKLdist:
                minKLdist=kld
                minKLshift=shift

        return (minKLdist,minKLshift)

    def maximumExpectedScore(self,other):
        "Compute the maximum expected score of self given other"

        self.computeFreqMatrix()
        other.computeFreqMatrix()


        maxEscore=-1e99
        for shift in range(len(other)):#range(1-len(self),len(other)):
            Es=self.maxExpectedScore(other,shift)
            #print shift,kld
            if Es>maxEscore:
                maxEscore=Es
                maxEshift=shift

        return (maxEscore,maxEshift)


    def __giveReverseComplementFreqs(self):
        "Return self.freq for reverse complement site"
        rc=[]
        for flist in self.freq:
            flist=flist[:]
            flist.reverse()
            rc.append(flist)
        rc.reverse()
        return rc

    def KLdistance(self,other,shift=0):
        "Computes the Kullback-Leibler distance min D(self||other) D(other||self) between two matrices."

        minKLdist=1e999
        for selfFreq in [self.freq,self.__giveReverseComplementFreqs()]:
            KLdistp=0.0
            KLdistq=0.0
            for (pList,qList,bgFreq) in zip(selfFreq,other.freq,[self.freqA,self.freqC,self.freqG,self.freqT]):
                if shift<0:
                    pList=[bgFreq]*(-shift)+pList
                else:
                    lenDelta=len(self)-len(other)
                    pList=pList[shift:]+[bgFreq]*(shift-lenDelta+1)
                qList=qList[:len(other)]
                pList=pList[:len(other)]
                #print "len",len(qList),len(pList)
                KLdistp+=reduce(operator.add,[ p*log2(p/q) for (p,q) in zip(pList,qList) ])
                KLdistq+=reduce(operator.add,[ q*log2(q/p) for (p,q) in zip(pList,qList) ])
            minKLdist=min(KLdistp,KLdistq,minKLdist)
            
        return minKLdist

    def maxExpectedScore(self,other,shift=0):
        "Computes maximum  (over other sites strand) expected score of self, given other."

        maxEscore=-1e999
        for otherFreq in [other.freq,other.__giveReverseComplementFreqs()]:
            Escore=0.0
            for (pList,scoreList,bgFreq) in zip(otherFreq,self.M_weight,[self.freqA,self.freqC,self.freqG,self.freqT]):
                if shift<0:
                    pList=[bgFreq]*(-shift)+pList
                else:
                    lenDelta=len(self)-len(other)
                    pList=pList[shift:]+[bgFreq]*(shift-lenDelta+1)
                scoreList=scoreList[:len(self)]
                pList=pList[:len(self)]
                #print "len",len(scoreList),len(pList)
                Escore+=reduce(operator.add,[ p*S for (p,S) in zip(pList,scoreList) ])

            maxEscore=max(Escore,maxEscore)

        #print self.name,other.name,maxEscore,shift
        return maxEscore


        

    def match(self,sequence):
        "matches matrix on sequence DOES NOT CURRENTLY WORK"
        #self.setBG(sequence)
        print "Match does not currently work"
        #return matrix.match(self.M_weight, sequence)

    def getTFBSbyAbsolute(self,sequence,cutoff):
        """Returns the hits that are better than cutoff"""
        #self.setBG(sequence)
        #self.initWeights()
        print "Max Score: %f Cutoff: %f"%(self.maxscore,cutoff)
        #seqIO=StringIO(sequence)
        #seqIO.seek(0)
        
        #self.initWeights(0)

        
        if self.maxscore>cutoff:
            #bg=matrix.BackGround(open("/home/kpalin/tyot/comparative/humanGenome/chr1.fa"))
            ret=_c_matrix.getTFBSwithBg(self.M_weight,sequence,cutoff,Matrix.backGround)
            #ret=matrix.getTFBS(self.M_weight,sequence,cutoff)
            if ret.has_key("NEXT_SEQ"):
                print "NEXT_SEQ",ret["NEXT_SEQ"]
                del ret["NEXT_SEQ"]
        else:
            ret={}
        return ret

    def getTFBSbyRatio(self, sequence, minscore_percent=0.1):
        """returns the hits, which are better then log2(minscore_percent*maxprod)
        
        These are possible Transcription Factor Binding Sites. The formula is equal
        to log2(minscore_percent)+maxscore"""
        minscore=log2(minscore_percent)
        #self.setBG(sequence)
        print "Max Score: %f Cutoff: %f"%(self.maxscore,minscore+self.maxscore)
        #seqIO=StringIO(sequence)
        #seqIO.seek(0)
        #seqIO=open("/home/kpalin/tyot/comparative/mabs/iso.tmp.fa")
        #seqIO.seek(1)
        ret=_c_matrix.getTFBSwithBg(self.M_weight, sequence, minscore+self.maxscore,Matrix.backGround)
        if not ret:
            ret={}
        if ret.has_key("NEXT_SEQ"):
            print "NEXT_SEQ",ret["NEXT_SEQ"]
            del ret["NEXT_SEQ"]
        return ret
    


def setPseudoCount(pcount=1.0):
    """Sets the amount of pseudocount"""
    if pcount<=0:
        print "Pseudocount must be non negative. Setting to one (1.0)!"
        pcount=1.0
    Matrix.pseudoCount=pcount

def setBGfreq(a,c,g,t):
    """Sets 0-order background.

    The parameters are frequences of a,c,g and t in the background sequence"""
    Matrix.backGround=None
    Matrix.freqA,Matrix.freqC,Matrix.freqG,Matrix.freqT=a,c,g,t


def setMarkovBackground(bg):
    """Set a markov Background.ds

    Markov background of k-order is represented as counts of (k+1)-grams
    in the background sequence"""
    print "Setting %dth order Markov background"%(bg.order)
    Matrix.backGround=bg

def getAllTFBS(sequence,cutoff,matlist,absoluteCutoff=None):
    "Get all TFBSs from one sequence."
    
    if not absoluteCutoff:
        cutoff=[log2(cutoff)+m.maxscore for m in matlist]

    if len(matlist)==0:
        ret={}
    else:
        Mat=matlist[:]
        ret=_c_matrix.getAllTFBSwithBg(Mat,sequence,cutoff,Matrix.backGround)


    return ret
