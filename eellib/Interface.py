"""Most of the Python functionality and gzip interface"""
# -*- coding: UTF-8 -*-

from Sequence import Sequences

from eellib import Matrix
import Output
import os,shutil
from tempfile import mktemp
import atexit
from glob import glob

import sys
if sys.platform!='win32':
    try:
        from gzip import GzipFile
    except ImportError:
        print "No gzip available."


#
# $Log$
# Revision 1.19  2004/12/22 11:14:24  kpalin
# Some fixes for better distributability
#
# Revision 1.18  2004/12/22 08:02:59  kpalin
# Hopefully more IO efficient TFBS search.
#
# Revision 1.17  2004/12/17 12:20:44  kpalin
# Changed the matrix dictionary key ordering
#
# Revision 1.16  2004/07/30 12:09:57  kpalin
# Working commands for multiple alignment.
#
# Revision 1.15  2004/04/14 07:48:07  kpalin
# Fixes to output and commands for MultiAlign
#
# Revision 1.14  2004/04/08 13:03:33  kpalin
# Updates on multiple alignment.
#
# Revision 1.13  2004/03/03 09:26:34  kpalin
# Added interface for multiple alignment.
#
# Revision 1.12  2004/02/23 12:23:52  kpalin
# Updates for per gene orthologous runs. Maybe litle multiple alignment.
#
# Revision 1.11  2004/02/13 09:40:19  kpalin
# Corrected bugs
#
# Revision 1.10  2004/02/11 09:38:13  kpalin
# Enabled memory saving features.
#
# Breaks 'more' command but saves a lot of memory.
#
# Revision 1.9  2004/02/05 10:31:47  kpalin
# Added extra garbage collection.
#
#


if 0:
    print "NOW IMPORTING ALIGN module"

from eellib import align

## try:
##     # Greedy multiple alignment written in python
##     from eellib import Multialign

##     # Exact multiD multiple alignment written in C
##     from eellib import multiAlign
## except ImportError:
##     pass

import sys,math


def memFormat(value):
    if value<2048:
        return "%dB"%(value)
    elif value<(1024*2048):
        return "%dkB"%(value/1024)
    elif value<(1024*1024*2048):
        return "%dMB"%(value/(1024*1024))
    #if value<(1024*1024*1024*2048):
    else:
        return "%dGB"%(value/(1024*1024*1024))

def timeFormat(value):
    s="%ds"%(value%60)
    value/=60
    if value>0:
        s="%dm %s"%(value%60,s)
        value/=60
    if value>0:
        s="%dh %s"%(value%24,s)
        value/=24
    if value>0:
        s="%dd %s"%(value,s)
    return s
    

def statReport():
    status={}
    try:
        statParts=open("/proc/self/stat").read().split()
        #print zip(range(len(statParts)),statParts)
        status["pid"]=int(statParts[0])
        status["comm"]=statParts[1]
        status["state"]=statParts[2]
        status["ppid"]=int(statParts[3])
        status["pgrp"]=int(statParts[4])
        status["session"]=int(statParts[5])
        status["tty_nr"]=int(statParts[6])
        status["tpgid"]=int(statParts[7])
        status["flags"]=long(statParts[8])
        status["minflt"]=long(statParts[9])
        status["majflt"]=long(statParts[10])
        status["cmajflt"]=long(statParts[11])
        status["utime"]=long(statParts[12])
        status["stime"]=long(statParts[13])
        status["cutime"]=long(statParts[14])
        status["cstime"]=long(statParts[15])
        status["priority"]=int(statParts[16])
        status["nice"]=int(statParts[17])  # number 18 is a placeholder
        status["itrealvalue"]=int(statParts[19])
        status["starttime"]=long(statParts[20])
        status["vsize"]=long(statParts[21])
        ## And a lot more. See: man proc
        return "Virtual memory size %s. Running time %s (%d seconds)."%(memFormat(status["vsize"]),timeFormat((status["cutime"]+status["cstime"])/100),(status["cutime"]+status["cstime"])/110.0)
    except Exception,e:
        pass

def memReport():
    try:
        import re
        statStr=open("/proc/self/status").read()
        m=re.search(r"VmSize:\s*(.+)",statStr)
        if m:
            print "Virtual memory size",m.group(1)
        else:
            print "Couldn't figure out status:\n",statStr
    except Exception:
        pass
def __sRep():
    print statReport()


atexit.register(memReport)
#atexit.register(__sRep)

try:
    from gc import collect
except ImportError:
    # Define a stubb if there is no Garbage Collection module
    def collect():
        pass

class Interface:
    """This class provides most of the functionality.

    Most of the UI commands are related to methods of this class.
    Especially align and getTFBS functions are important."""
    def __init__(self):
        self.outputted=0
        self.resetMatrices(self)
        self.resetSequences(self)
        self.__comp={}
        self.alignment=None


    def multiAlignGreedy(self,arglist):
        """Arguments: pairwiseGFFfiles
        Join the given pairwise alignment GFF files to multiple alignment."""
        #self.malignment=Multialign.MultipleAlignment()

        for fileGlob in arglist:
            filenames=glob(fileGlob)
            if len(filenames)==0:
                print "Can't find",fileGlob
            for fileName in filenames:
                print "Reading",fileName
                self.malignment.addGFFfile(fileName)
        print "All %d files added. Doing the alignment"%(len(filenames)),filenames
        self.malignment.multiAlign()



    def multiAlign(self,arglist):
        """Arguments: [filename[,num_of_align,[lambda[,xi[,mu[,nu,[,nuc_per_rotation]]]]]]]
Computes multiple alignment of the BS or optionally the BS from a gff file
filename specifies a file in gff format is you want to be aligned
num_of_align        specifies how many alignments you want. (Default 3)
lambda   Bonus factor for hit (Default 2)
xi       Penalty factor for rotation (Default 1.0)
mu       Penalty factor for average distance between sites (Default 0.5)
nu       Penalty factor for distance difference between sites (Default 1.0)
nuc_per_rotation    specifies how many nucletides there are per rotation. (Default 10.4)
If you want to skip a argument just  write '.' for it.
If you use '.' as filename the local data are aligned."""
        try:
            [filename, num_of_align, Lambda, xi,
             mu, nu, nuc_per_rotation]=arglist + ['.']*(7-len(arglist))
            
            if num_of_align=='.':
                num_of_align=3
            if Lambda=='.':
                Lambda=2.0
            if xi=='.':
                xi=1.0
            if mu=='.':
                mu=0.5
            if nu=='.':
                nu=1.0
            if nuc_per_rotation=='.':
                nuc_per_rotation=10.4


            data=[x.split() for x in open(filename).readlines() ]
            self.alignment=multiAlign.MultiAlignment(data,int(num_of_align),
                                                     float(Lambda), float(xi),
                                                     float(mu), float(nu),float(nuc_per_rotation))
            if not self.alignment:
                print "No multiple alignment for a reason or an other"
            else:
                print "Done"
#                print "goodAlign=",self.malignment.nextBest()
        except ValueError:
            print "Error: unallowed arguments passed to 'multipleAlign'"

        
    def showMultiAlign(self,arglist):
        """Arguments: [minPairs]
        Outputs the multiple alignment to standard output.
        If an integer minPairs is given, no sites aligned with less than that
        number of pairwise alignments, is reported."""
        if not hasattr(self,"malignment"):
            return
        if len(arglist)>0:
            minPairs=int(arglist[0])
        else:
            minPairs=0
        #for i in [self.malignment[0]]:
        for i in self.malignment:
            i.setAlnLimit(minPairs)
            #if not i or len(i)<2:continue
            if len([x for x in i.seqs if x in self.seq.getNames()])==len(i.seqs):
                i.strAln(self.seq)
            print str(i)
            print "\n"


    def saveMultiAlign(self,arglist):
        """Arguments: filename [minPairs]
        Outputs the multiple alignment to file 'filename'.
        If an integer minPairs is given, no sites aligned with less than that
        number of pairwise alignments, is reported."""
        fname=arglist[0]
        try:
            minPairs=int(arglist[1])
        except (IndexError,ValueError):
            minPairs=0
            
        if not hasattr(self,"malignment"):
            print "No multiple alignment to save!"
            return
        m="w"
        for i in self.malignment:
            i.setAlnLimit(minPairs)
            #if len(i)<2:continue
            if len([x for x in i.seqs if x in self.seq.getNames()])==len(i.seqs):
                i.strAln(self.seq)
            Output.savealign(str(i)+"\n",fname,m)
            m="a"

    def resetMatrices(self):
        "resets the list of matrices"
        self.matlist=[]
        collect()

    def addMatrix(self, filenames):
        "adds new matrices from given files"
        for f in filenames:
            try:
                m=Matrix.Matrix(f)
                print "adding %s, Info=%2.4g:"%(f,m.InfoContent)
                self.matlist.append(m)
            except ValueError:
                print "could not read",f
            except IOError, (errno, strerror):
                print "%s: %s" % (strerror, f)


            
    def printMatrices(self):
        "prints matrices to standard out"
        for i in range(len(self.matlist)):
            print "Matrix No %d (info=%g) %s"%(i,self.matlist[i].InfoContent,self.matlist[i].name)
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
        collect()

    def addSequence(self, filenames):
        "adds sequences from files in FASTA format"
        for f in filenames:
            self.seq.addSequence(f)
            
    def removeSequence(self, name):
        "removes sequences given by name"
        self.seq.removeSequence(name)
        collect()
    
    def getSeqNames(self):
        "returns the sequence names"
        return self.seq.getNames()


    def getTFBSAbsolute(self,cutoff=9.0):
        "computes the possible TFBS"
        Interface.getTFBS(self,cutoff,1)
 
    def getTFBS(self, bound=0.1 , absCutoff=None):
        "computes the possible TFBS"
        if hasattr(self,"tempFileName"):
            os.remove(self.tempFileName)
            del(self.tempFileName)

        self.__comp={}
        duration = len(self.matlist) * len(self.seq.getNames())
        progress= 0.0
        totalMatches=0
        seqnumber=0
        for name in self.seq.getNames():
            seqnumber +=1
            print "Matching Sequence.",seqnumber,"of",len(self.seq.getNames())
            self.__comp[name]={}
            self.__comp[name]=Matrix.getAllTFBS(self.seq.sequence(name),
                                                  bound,self.matlist,absCutoff)
#                self.__comp[name][m]=m.getTFBSbyRatio(self.seq.sequence(name),
#                                                      bound)
            try:
                foundMatches=sum([len(x) for x in self.__comp[name].values()])
                totalMatches+=foundMatches
                print "Found %d matches\n"%(foundMatches)
            except TypeError:
                print "name=",name
                print "self.__comp=",self.__comp
                #self.__gff=Output.get(self.__comp).split('\n')
            if totalMatches>50000:
                self.storeTmpGFF()
                self.__comp[name]={}
                totalMatches=0
        if hasattr(self,"tempFileName"):
            self.finalTmpGFF()

    def storeTmpGFF(self):
        "Open temporary file and store BS data to it"
        if not hasattr(self,"tempFile"):
            self.tempFileName=mktemp(".gff.gz")
            try:
                self.tempFile=GzipFile(self.tempFileName,"w")
            except NameError:
                self.tempFileName=self.tempFileName[:-3]
                self.tempFile=open(self.tempFileName,"w")
            def condRemoveTmp(tmpName):
                try:
                    os.remove(tmpName)
                except OSError:
                    pass
            atexit.register(condRemoveTmp,self.tempFileName)
            print "Storing temporary file",self.tempFileName
        print "Writing temporary file"
        outData=Output.get(self.__comp)
        if len(outData)>0:
            self.tempFile.write(outData)
            self.tempFile.flush()
            self.__comp={}

    def finalTmpGFF(self):
        "Store BS data to temporary file and close it"
        self.storeTmpGFF()
        self.tempFile.close()
        del(self.tempFile)
        collect()
        
            
    def savematch(self, filename=''):
        "Saves the results. (possibly gziped)"
        # To view the gff file numbered in the same order, use:
        #grep  ENSG tmp.gff |sort --key=5n,5 --key=6n |nl -v0>tmp.sort.ensg.gff
        if hasattr(self,"tempFileName"):
            print "Storing file in gziped format."
            if not filename[-2:]=="gz":
                filename=filename+".gz"
            try:
                shutil.copy(self.tempFileName,filename)
                os.remove(self.tempFileName)
                del(self.tempFileName)
            except IOError,e:
                print e
                filename=self.tempFileName
            return filename
        else:
            return Output.savematch(self.__comp, filename)

    def showmatch(self):
        "Prints the results to standard out"
        if hasattr(self,"tempFileName"):
            try:
                tempFile=GzipFile(self.tempFileName,"r")
            except NameError:
                tempFile=open(self.tempFileName,"r")
            print tempFile.read()
            tempFile.close()
        else:
            Output.showmatch(self.__comp)

    def showMoreAlignments(self,count=1):
        """Gets more alignments but doesn't display them"""
        self.moreAlignments(count)
        #print Output.formatalign(self.alignment,self.seq),
        

    def quit(self):
        "Exits the program"
        print "Exiting the program"
        sys.exit()


    def moreAlignments(self,num_of_align=1):
        """Fetch more alignments from previously run alignment matrix"""
        if not self.alignment or not hasattr(self.alignment,"memSaveUsed"):
            print self.alignment
            return
        for i in range(num_of_align):
            if self.alignment.memSaveUsed==1 and self.alignment.askedResults<=len(self.alignment.bestAlignments):
                print "Can't give more alignments. Don't remember those"
                break
            goodAlign=self.alignment.nextBest()
            if not goodAlign:
                break
            goodAlign.reverse() # For old times sake
            #               0 1 2     3      4[0]   4[1]   5[0]   5[1]  6
            # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
            self.alignment.bestAlignments.append(goodAlign)
            
        
    def align(self, filename='.', num_of_align=3,
              Lambda=1.0, xi=1.0, mu=0.5,nu=1.0, nuc_per_rotation=10.4):
        "aligns the computed BS or some given in file"
        if hasattr(self,"alignment"):
            del(self.alignment)
            
        if filename=='.' and  hasattr(self,"tempFileName"):
            filename=self.tempFileName

        collect()
        if filename=='.':
            if len(self.__comp)==0:
                return "No binding sites"
            else:
                self.alignment=align.aligndata(Output.get(self.__comp),num_of_align,
                                               Lambda, xi,
                                               mu, nu,nuc_per_rotation)

        else:
            self.alignment= align.alignfile(filename, num_of_align,Lambda,
                                            xi, mu, nu,nuc_per_rotation)

        self.moreAlignments(num_of_align)
        if self.alignment:
            #Interface.showalignSTDO(self)
            assert self.alignment.x_name!=self.alignment.y_name
            print "Used time %g sec."%(self.alignment.secs_to_align)
            return 1
        else:
            return None

    def savealignGFF(self,filename=""):
        "Saves the results in GFF format"
        return Output.savealign(Output.formatalignGFF(self.alignment), filename)

    def savealignAnchor(self,filename=""):
        "Saves the results in Anchor format"
        return Output.savealign(Output.formatalignCHAOS(self.alignment), filename)
        
    def savealign(self, filename=''):
        "Saves the results"
        return Output.savealign(Output.formatalign(self.alignment,self.seq), filename)

    def showalignSTDO(self):
        "Prints the alignment to standart out"
        print Output.formatalign(self.alignment,self.seq),

    def about(self):
        "Information about the program"
        return """Enhancer Element Locator EEL
Version $Version$
Contact: Kimmo Palin (kimmo.palin@helsinki.fi)

Big thanks to Matthias Berg (bergm@web.de) for coding.
"""
#         Uni 18, Zi. 2220
#         D-66123 Saarbruecken
#         Germany


