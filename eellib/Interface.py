from Matrix import Matrix
from Sequence import Sequences
import Output
import os,shutil
from tempfile import mktemp
try:
    from gzip import GzipFile
except ImportError:
    print "No gzip available."

if 0:
    print "NOW IMPORTING ALIGN module"
import align
import sys,math

class Interface:
    #This class provides the functionality which is used by user interfaces
    def __init__(self):
        self.outputted=0
        self.resetMatrices(self)
        self.resetSequences(self)
        self.__comp={}
        self.alignment=None

    def resetMatrices(self):
        "resets the list of matrices"
        self.matlist=[]

    def addMatrix(self, filenames):
        "adds new matrices from given files"
        for f in filenames:
            try:
                m=Matrix(f)
                print "adding",f
                self.matlist.append(m)
            except ValueError:
                print "could not read",f
            except IOError, (errno, strerror):
                print "%s: %s" % (strerror, f)


            
    def printMatrices(self):
        "prints matrices to standard out"
        for i in range(len(self.matlist)):
            print "Matrix No %d  %s"%(i,self.matlist[i].name)
            self.matlist[i].draw()

    def printMatrixWeights(self):
        "prints matrix weights to standard out"
        for i in range(len(self.matlist)):
            try:
                bgStr="Order %d Markov"%(self.matlist[i].backGround.order)
            except AttributeError:
                bgStr="(%0.2f,%0.2f,%0.2f,%0.2f)"%(self.matlist[i].freqA,self.matlist[i].freqC,self.matlist[i].freqG,self.matlist[i].freqT)
            print "Matrix No %d %s bg=%s pCount=%g maxscore=%0.2f"%(i,self.matlist[i].name,bgStr,self.matlist[i].pseudoCount,self.matlist[i].maxscore)
            self.matlist[i].drawWeights()

    def removeMatrix(self, index):
        "removes matrix given by index"
        self.matlist.pop(index)

    def resetSequences(self):
        "resets the sequences"
        self.seq=Sequences()

    def addSequence(self, filenames):
        "adds sequences from files in FASTA format"
        for f in filenames:
            self.seq.addSequence(f)
            
    def removeSequence(self, name):
        "removes sequences given by name"
        self.seq.removeSequence(name)
    
    def getSeqNames(self):
        "returns the sequence names"
        return self.seq.getNames()


    def getTFBSAbsolute(self,cutoff=9.0):
        "computes the possible TFBS"
        Interface.getTFBS(self,cutoff,1)
 
    def getTFBS(self, bound=0.1 , absCutoff=None):
        "computes the possible TFBS"
        if hasattr(self,"tempFile"):
            del(self.tempFile)
            os.remove(self.tempFileName)
            del(self.tempFileName)

        self.__comp={}
        duration = len(self.matlist) * len(self.seq.getNames())
        progress= 0.0
        totalMatches=0
        matrixnumber=0
        for m in self.matlist:
            matrixnumber +=1
            print "Matching matrix no.",matrixnumber,"of",len(self.matlist)
            self.__comp[m]={}
            for name in self.seq.getNames():
                print "Matching matrix %s on sequence %s"%(m.getName(),name)
                print "Progress: %2.2f %%" % (100*progress/duration)
                progress +=1.0
                if absCutoff:
                    self.__comp[m][name]=m.getTFBSbyAbsolute(self.seq.sequence(name),
                                                             bound)
                else:
                    self.__comp[m][name]=m.getTFBSbyRatio(self.seq.sequence(name),
                                                          bound)
                try:
                    totalMatches+=len(self.__comp[m][name])
                    print "Found %d matches\n"%(len(self.__comp[m][name]))
                except TypeError:
                    print "m=",m
                    print "name=",name
                    print "self.__comp=",self.__comp
                    #self.__gff=Output.get(self.__comp).split('\n')
                if hasattr(self,"tempFile") or totalMatches>50000:
                    self.storeTmpGFF()
        if hasattr(self,"tempFile"):
            self.storeTmpGFF()

    def storeTmpGFF(self):
        if not hasattr(self,"tempFile"):
            self.tempFileName=mktemp(".gff.gz")
            try:
                self.tempFile=GzipFile(self.tempFileName,"w")
            except NameError:
                self.tempFile=self.tempFileName[:-3]
                self.tempFile=open(self.tempFileName,"w")
            print "Storing gziped temporary file",self.tempFileName
        self.tempFile.write(Output.get(self.__comp))
        self.tempFile.flush()
        self.__gff=""
        
            
    def savematch(self, filename=''):
        "Saves the results."
        # To view the gff file numbered in the same order, use:
        #grep  ENSG tmp.gff |sort --key=5n,5 --key=6n |nl -v0>tmp.sort.ensg.gff
        if hasattr(self,"tempFile"):
            print "Storing file in gziped format."
            self.tempFile.close()
            if not filename[-2:]=="gz":
                filename=filename+".gz"
            try:
                shutil.copy(self.tempFileName,filename)
                os.remove(self.tempFileName)
                del(self.tempFile)
                del(self.tempFileName)
            except IOError,e:
                print e
                filename=self.tempFileName
            return filename
        else:
            return Output.savematch(self.__comp, filename)

    def showmatch(self):
        "Prints the results to standard out"
        if hasattr(self,"tempFile"):
            self.tempFile.seek(0)
            print self.tempFile.read()
        else:
            Output.showmatch(self.__comp)

    def showMoreAlignments(self,count=1):
        self.moreAlignments(count)
        #print Output.formatalign(self.alignment,self.seq),
        

    def quit(self):
        "Exits the program"
        sys.exit()


    def moreAlignments(self,num_of_align=1):
        if not self.alignment:
            return
        for i in range(num_of_align):
            goodAlign=self.alignment.nextBest()
            if not goodAlign:
                break
            goodAlign.reverse() # For old times sake

            # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
            self.alignment.bestAlignments.append(goodAlign)
            
        
    def align(self, filename='.', num_of_align=3,
              Lambda=1.0, xi=1.0, mu=0.5,nu=1.0, nuc_per_rotation=10.4):
        "aligns the computed BS or some given in file"
        if hasattr(self,"alignment"):
            del(self.alignment)
        if filename=='.' and  hasattr(self,"tempFile"):
            filename=self.tempFileName
        elif filename=='.' and len(self.__comp)==0:
            return "No binding sites"
        elif filename==".":
            self.alignment=align.aligndata(Output.get(self.__comp),
                                             num_of_align, Lambda, xi,
                                             mu, nu,nuc_per_rotation)
        else:
            self.alignment= align.alignfile(filename, num_of_align, Lambda,
                                              xi, mu, nu,nuc_per_rotation)

        self.moreAlignments(num_of_align)
        if self.alignment:
            Interface.showalignSTDO(self)
            assert self.alignment.x_name!=self.alignment.y_name
            return 1
        else:
            return None

    def savealignGFF(self,filename=""):
        "Saves the results in GFF format"
        return Output.savealign(Output.formatalignGFF(self.alignment), filename)
        
    def savealign(self, filename=''):
        "Saves the results"
        return Output.savealign(Output.formatalign(self.alignment,self.seq), filename)

    def showalignSTDO(self):
        "Prints the alignment to standart out"
        print Output.formatalign(self.alignment,self.seq),

    def about(self):
        "Information about the program"
        return """Matthias Berg's alignment of binding sites
Version 1.1
Contact: Matthias Berg
         Uni 18, Zi. 2220
         D-66123 Saarbruecken
         Germany

         email: bergm@web.de"""

