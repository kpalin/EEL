"""Python interface and logic to TFBS matrix matching"""
import operator
import string
import matrix
import math


#
# $Log$
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
    def __init__(self,filename):
        "reads matrix from file"
        self.pseudoCount=1.0

        self.backGround=None
        File=open(filename,'r')
        self.name=filename
        self.LLMatrix=filter(lambda x:len(x),[[string.atoi(entry)
                        for entry in line.split()]
                       for line in File.read().split('\n')])
        File.close()
        self.setBGfreq(0.25,0.25,0.25,0.25)
        #self.draw()


    def setMarkovBackground(self,bg):
        """Set a markov Background.

        Markov background of k-order is represented as counts of (k+1)-grams
        in the background sequence"""
        self.backGround=bg
        self.initWeights()

    def setBGfreq(self,a,c,g,t):
        """Sets 0-order background.

        The parameters are frequences of a,c,g and t in the background sequence"""
        self.backGround=None
        self.freqA,self.freqC,self.freqG,self.freqT=a,c,g,t
        self.initWeights()

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

    def setPseudoCount(self,pcount=1.0):
        """Sets the amount of pseudocount"""
        if pcount<=0:
            print "Pseudocount must be non negative. Setting to one (1.0)!"
            pcount=1.0
        self.pseudoCount=pcount
        self.initWeights()

    def initWeights(self):
        """Helper to initialize the matrix weights for 0- or higher order background models"""
        self.M_weight=[]
        if self.backGround:
            self.positiveWeights()
        else:
            self.trivialWeights()

    def positiveWeights(self):
        """Update weights for higher order markov background.

        Update weights only for positive probability.
        Background is taken care elsewhere"""

        sum=reduce(lambda s,row:map(operator.add,s,row),self.LLMatrix,[self.pseudoCount]*len(self.LLMatrix[0]))


        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqA*self.pseudoCount))/tot),
                            self.LLMatrix[0],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*self.pseudoCount))/tot),
                            self.LLMatrix[1],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*self.pseudoCount))/tot),
                            self.LLMatrix[2],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqT*self.pseudoCount))/tot),
                            self.LLMatrix[3],sum))

        # Don't know how to compute maxscore with markov background. It varies.
        #maxscore is the highest reachable score (disregarding BG)
        #self.maxscore=0.0
        #for i in range(len(self.M_weight[0])):
        #    self.maxscore-=max(self.M_weight[0][i],
        #                         self.M_weight[1][i],
        #                         self.M_weight[2][i],
        #                         self.M_weight[3][i])

        

    def trivialWeights(self):
        """Update weights according to 0-order background"""

        # Sum of the columns.
        sum=reduce(lambda s,row:map(operator.add,s,row),self.LLMatrix,[self.pseudoCount]*len(self.LLMatrix[0]))


        # Frequencies of a nucleotide
        freq=[ [ (x+bgFreq*self.pseudoCount)/(tot) for x,tot in zip(matrix,sum) ] \
               for matrix,bgFreq in zip(self.LLMatrix,(self.freqA,self.freqC,self.freqG,self.freqT)) ]

        # Compute information content of the motif
        InfoContent=[[ f*log2(f/bgF) for f in fLine] for (fLine,bgF) in zip(freq,(self.freqA,self.freqC,self.freqG,self.freqT)) ]


        self.InfoContent=reduce(operator.add,reduce(lambda x,y:x+y,InfoContent,[]))
        # Compute weights.
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqA*self.pseudoCount))/
                                               (tot*self.freqA)),
                            self.LLMatrix[0],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqC*self.pseudoCount))/
                                               (tot*self.freqC)),
                            self.LLMatrix[1],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqG*self.pseudoCount))/
                                               (tot*self.freqG)),
                            self.LLMatrix[2],sum))
        self.M_weight.append(map(lambda x,tot: log2((x+(self.freqT*self.pseudoCount))/
                                               (tot*self.freqT)),
                            self.LLMatrix[3],sum))


        #maxscore is the highest reachable score (which respect to the BG)
        self.maxscore=0.0
        for i in range(len(self.M_weight[0])):
            self.maxscore+=max(self.M_weight[0][i],
                                 self.M_weight[1][i],
                                 self.M_weight[2][i],
                                 self.M_weight[3][i])



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
            ret=matrix.getTFBSwithBg(self.M_weight,sequence,cutoff,self.backGround)
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
        ret=matrix.getTFBSwithBg(self.M_weight, sequence, minscore+self.maxscore,self.backGround)
        if not ret:
            ret={}
        if ret.has_key("NEXT_SEQ"):
            print "NEXT_SEQ",ret["NEXT_SEQ"]
            del ret["NEXT_SEQ"]
        return ret
    



def getAllTFBS(sequence,cutoff,matlist,absOrRat=None):
    "Get all TFBSs from one sequence."
    
    if absOrRat:
        matlist=[x for x in matlist if x.maxscore>cutoff]

    else:
        cutoff=[log2(cutoff)+m.maxscore for m in matlist]

    if len(matlist)==0:
        ret={}
    else:
        Mat=matlist[:]
        BG=matlist[0].backGround
        ret=matrix.getAllTFBSwithBg(Mat,sequence,cutoff,BG)


    return ret
