"""Most of the Python functionality and gzip interface"""
# -*- coding: UTF-8 -*-

from Sequence import Sequences

from eellib import Matrix

import Output
import os,shutil
from tempfile import mktemp
import atexit
from glob import glob
from time import time
import operator

class CommandError(Exception):
    "Error that occured during processing of a command"
    def __init__(self,message,parent=None):
        """Initialize the error.

        message is the message output after raising the error.
        parent is the original exception that resulted in the error."""
        self.message=message
        self.parent=parent

    def __str__(self):
        return self.message

import sys
if sys.platform!='win32':
    try:
        from gzip import GzipFile
    except ImportError:
        print "No gzip available."


#
# $Log$
# Revision 1.46  2008/04/01 05:47:20  kpalin
# Fixed an indentation.
#
# Revision 1.45  2008/02/29 08:53:32  kpalin
# Should provide output even when the memory save was used
#
# Revision 1.44  2007/12/13 09:33:31  kpalin
# Fixed PFM file name pretty printing and increased temporary GFF file
# creation limit.
#
# Revision 1.43  2007/08/09 05:59:59  kpalin
# Fixed setBGFreq to handle missing sample sequence request gracefully.
#
# Revision 1.42  2007/08/09 05:58:02  kpalin
# Fixed setBGFreq to use sample sequence
#
# Revision 1.41  2007/06/05 12:40:00  kpalin
# Dont' remember
#
# Revision 1.40  2007/02/08 09:21:34  kpalin
# Added possibility of setting simple background from example.
#
# Revision 1.39  2006/12/08 09:51:34  kpalin
# Added E-values and TFBS p-values
#
# Revision 1.38  2006/11/13 13:03:22  kpalin
# Added E-model RMSE
#
# Revision 1.37  2006/11/13 12:34:48  kpalin
# Added TFBS search p-value cutoff and Escore computation.
#
# Revision 1.36  2006/10/26 08:52:37  kpalin
# Workaround for old python v2.2
#
# Revision 1.35  2006/10/18 07:23:44  kpalin
# Added suboptimalsDownTO
#
# Revision 1.34  2006/08/31 11:19:22  kpalin
# Align pre-given sequences and get TFBS only if needed.
#
# Revision 1.33  2006/08/31 10:08:36  kpalin
# More Debug output and TF distance calculations.
#
# Revision 1.32  2006/08/14 09:57:30  kpalin
# Added randomization function and cleaned diagnostic outputs.
#
# Revision 1.31  2006/04/05 08:30:07  kpalin
# Regular expressions for sequence removal and commands for multiple alignment.
#
# Revision 1.30  2005/06/10 06:55:08  kpalin
# TIming for matrix matching.
#
# Revision 1.29  2005/05/19 07:49:35  kpalin
# Merged Waterman-Eggert style suboptimal alignments and
# SNP matching.
#
# Revision 1.28.2.3  2005/05/09 07:14:07  kpalin
# Catch exception with faulty savealignGFF
#
# Revision 1.28.2.2  2005/04/12 09:10:44  kpalin
# Suboptimals is the default way of getting new results.
#
# Revision 1.28.2.1  2005/03/31 13:30:54  kpalin
# Added command 'suboptimal' which is like 'more' but gives
# real suboptimal results instead of next best from the alignment matrix.
#
# Revision 1.28  2005/03/22 13:17:13  kpalin
# Merged some fixes surfacing from testing the public version.
#
# Revision 1.22.2.3  2005/03/22 12:19:26  kpalin
# Fixed the problems that came up with testing.
#
# Revision 1.22.2.2  2005/03/08 11:14:03  kpalin
# Fix matrix/background position dependency issue.
#
# Revision 1.26  2005/03/08 10:51:34  kpalin
# Merging matrix/background change to distribution version.
#
# Revision 1.25  2005/03/08 10:38:23  kpalin
# Fixed markov background setting to keep effect also after adding
# more matricies.
#
# Revision 1.24  2005/02/24 11:36:53  kpalin
# Added handling of the site annotations.
#
# Revision 1.23  2005/02/21 09:50:41  kpalin
# Fixed a bug conserning beatifying matrix names with similar filenames.
#
# Revision 1.22  2005/01/13 13:16:42  kpalin
# Moved the requesting of sequences to be aligned to Python side
# of stuff. Much better.
#
# Revision 1.21  2005/01/12 13:34:55  kpalin
# Added Tkinter/Tix Graphical user interface and command -no-gui to
# avoid it.
#
# Revision 1.20  2005/01/07 13:41:25  kpalin
# Works with py2exe. (windows executables)
#
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


try: #Workaround for python 2.2
    sum([1,2])
except NameError:
    import operator
    sum=lambda x:reduce(operator.add,x)

if 0:
    print "NOW IMPORTING ALIGN module"

from eellib import align

try:
    # Greedy multiple alignment written in python
    from eellib import Multialign
    
    # Exact multiD multiple alignment written in C
    from eellib import multiAlign
except ImportError:
    pass

try:
    from eellib import shortMultiAlign
    from eellib import multiFromPairwise
    from eellib import treeMultiAlign
except ImportError:
    pass

from eellib import _c_matrix

import sys,math
import os.path

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
    value=int(value)
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
    "Return some status information about the program."
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
    "Report memory consumption"
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
        self.resetMatrices()
        self.resetSequences()
        self.resetTFBS()
        self.alignment=None


    def resetTFBS(self):
        self.__comp_bound=-1.0
        self.__comp_absCutoff=None
        self.__comp={}
        
    def showFileList(self,files):
        "Display the list of files on display device"
        for fileName in files:
            self.show(fileName)

    def show(self,*text):
        "Display the parameters on display device (stdout)"
        sys.stdout.writelines(" ".join(map(str,text+('\n',))))
            
    def head(self,sitesCount):
        "Keep only the first sitesCount binding sites for each sequence."

        allSites={}
        allLimits={}
        for sequence in self.__comp.keys():
            allSites[sequence]=[]
            for sites in self.__comp[sequence].values():
                allSites[sequence].extend([x[0] for x in sites])
            allSites[sequence].sort()
            if allSites[sequence]>sitesCount:
                allLimits[sequence]=allSites[sequence][sitesCount]+1
            else:
                allLimits[sequence]=allSites[-1]+1
            for matrix in self.__comp[sequence].keys():
                self.__comp[sequence][matrix]=dict([(x,y) for (x,y) in self.__comp[sequence][matrix].items() if x[0]<allLimits[sequence]])
        self.show("Cut site sequences to length %d\n"%(sitesCount))
        
    def multiAlignGreedy(self,arglist):
        """Arguments: pairwiseGFFfiles
        Join the given pairwise alignment GFF files to multiple alignment."""
        #self.malignment=Multialign.MultipleAlignment()

        for fileGlob in arglist:
            filenames=glob(fileGlob)
            if len(filenames)==0:
                self.show( "Can't find",fileGlob)
            for fileName in filenames:
                self.show("Reading",fileName)
                self.malignment.addGFFfile(fileName)
        self.show("All %d files added. Doing the alignment"%(len(filenames)),filenames)
        self.malignment.multiAlign()

    def multiFromPairwise(self,arglist):
        """Arguments: filename,pwfilename[,gaplimit[,numofalign[,buffer[,lambda[,xi[,mu[,nu,[,nuc_per_rotation]]]]]]
Computes multiple alignments from places selected from pairwise alignments of the BS from a gff file
filename   specifies a file in gff format is you want to be aligned
pwfilename specifies a file in alignment gff format where the pairwise alignments are
gaplimit   specifies how many sequences can have gaps at the same point (Default 0)
numofalign        specifies how many alignments you want. (Default 1)
buffer     Specifies the size of buffer to be added to the found areas before multiple aligning in nukleotides (Default 1000)
lambda     Bonus factor for hit (Default 2)
xi         Penalty factor for rotation (Default 1.0)
mu         Penalty factor for average distance between sites (Default 0.5)
nu         Penalty factor for distance difference between sites (Default 1.0)
nuc_per_rotation    Specifies how many nucletides there are per rotation. (Default 10.4)
If you want to skip a argument just  write '.' for it."""
#        self.show("Aloitetaan shortMultiAlign")
        try:
            [filename, pwfilename, gaplimit, numofalign, nbuffer, Lambda, xi,
             mu, nu, nuc_per_rotation]=arglist + ['.']*(10-len(arglist))
            
            if gaplimit=='.':
                gaplimit=0
            if numofalign=='.':
                numofalign=1
            if nbuffer=='.':
                nbuffer=1000
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


            data=[x.strip().split() for x in open(filename).readlines() ]
            while(len(data[-1])==0):
                del data[-1]
            
	    pairwise=[x.strip().split() for x in open(pwfilename).readlines() ]
            while(len(pairwise[-1])==0):
                del pairwise[-1]

            self.alignment=multiFromPairwise.MultiFromPairwise(data,pairwise,int(gaplimit),int(numofalign),
                                                     int(nbuffer), float(Lambda), float(xi),
                                                     float(mu), float(nu),float(nuc_per_rotation))
            if not self.alignment:
                self.show("No multiple alignment for a reason or an other")
            else:
                for alnNo,goodAlign in enumerate(self.alignment.bestAlignments):
                    goodAlign.reverse();
                self.show("Done")
                #self.show(str(len(self.alignment.bestAlignments)))
                self.show("Used time %g sec."%(self.alignment.secs_to_align))
        
	except ValueError,e:
            raise CommandError("Error:",e)
            print "Error: unallowed arguments passed to 'multiFromPairwise'"


    def shortMultiAlign(self,arglist):
        """Arguments: filename[,gaplimit[,numofalign[,lambda[,xi[,mu[,nu,[,nuc_per_rotation]]]]]]
Computes multiple alignment of the BS from a gff file
filename specifies a file in gff format is you want to be aligned
gaplimit specifies how many sequences can have gaps at the same point (Default 0)
numofalign        specifies how many alignments you want. (Default 3)
lambda   Bonus factor for hit (Default 2)
xi       Penalty factor for rotation (Default 1.0)
mu       Penalty factor for average distance between sites (Default 0.5)
nu       Penalty factor for distance difference between sites (Default 1.0)
nuc_per_rotation    Specifies how many nucletides there are per rotation. (Default 10.4)
If you want to skip a argument just  write '.' for it."""
#        self.show("Aloitetaan shortMultiAlign")
        try:
            [filename, gaplimit, numofalign, Lambda, xi,
             mu, nu, nuc_per_rotation]=arglist + ['.']*(8-len(arglist))
            
            if gaplimit=='.':
                gaplimit=0
            if numofalign=='.':
                numofalign=3
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


            data=[x.strip().split() for x in open(filename).readlines() ]
            while(len(data[-1])==0):
                del data[-1]

            self.alignment=shortMultiAlign.ShortMultiAlignment(data,int(gaplimit),int(numofalign),
                                                     float(Lambda), float(xi),
                                                     float(mu), float(nu),float(nuc_per_rotation))
            if not self.alignment:
                self.show("No multiple alignment for a reason or an other")
            else:
                for alnNo,goodAlign in enumerate(self.alignment.bestAlignments):
                    goodAlign.reverse()
                self.show("Done")
                #self.show(str(len(self.alignment.bestAlignments)))
                self.show("Used time %g sec."%(self.alignment.secs_to_align))
        
	except ValueError,e:
            raise CommandError("Error:",e)
            print "Error: unallowed arguments passed to 'shortMultipleAlign'"

    def treeMultiAlign(self,arglist):
        """Arguments: filename,treefilename[,method[,k[,numofalign[,lambda[,xi[,mu[,nu,[,nuc_per_rotation[,pwfiles]]]]]]]]]
Computes multiple alignment of the BS from a gff file
filename specifies a file in gff format is you want to be aligned
treefilename specifies a file in Newick tree format with the phylogenetic tree
method specifies which method of scoring is used. 1=Sum of pairs, 2=Score generalized from pairwise, 3=Score based on representative sequences. (Default 1)
k specifies how many local alignments are kept at each level (Default 100)
numofalign        specifies how many alignments you want. (Default 3)
lambda   Bonus factor for hit (Default 2)
xi       Penalty factor for rotation (Default 1.0)
mu       Penalty factor for average distance between sites (Default 0.5)
nu       Penalty factor for distance difference between sites (Default 1.0)
nuc_per_rotation    Specifies how many nucletides there are per rotation. (Default 10.4)
pwfiles	 If the needed pairwise alignments are calculated beforehand, the path to the directory they are in.
If you want to skip a argument just  write '.' for it."""
        try:
            [filename, treefilename, method, k, numofalign, Lambda, xi,
             mu, nu, nuc_per_rotation, pwfiles]=arglist + ['.']*(11-len(arglist))

            if k=='.':
                k=100
            if numofalign=='.':
                numofalign=3
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
            if method=='.':
                method=1
            if pwfiles=='.':
                pwfiles=""

            data=[x.strip().split() for x in open(filename).readlines() ]
            while(len(data[-1])==0):
                del data[-1]

#            tree=[x.strip().split() for x in open(treefilename).readlines() ]
            x=open(treefilename).readline()
            x.strip()
            tree=x

            self.alignment=treeMultiAlign.TreeMultiAlignment(data,tree,int(k),int(numofalign),int(method),
                                                     float(Lambda), float(xi),
                                                     float(mu), float(nu),float(nuc_per_rotation),pwfiles)
            if not self.alignment:
                self.show("No multiple alignment for a reason or an other")
            else:
#                for alnNo,goodAlign in enumerate(self.alignment.bestAlignments):
#                    goodAlign.reverse()
                self.show("Done")
                #self.show(str(len(self.alignment.bestAlignments)))
                self.show("Used time %g sec."%(self.alignment.secs_to_align))
        
	except ValueError,e:
            raise CommandError("Error:",e)
            print "Error: unallowed arguments passed to 'treeMultipleAlign'"


    def nodeAlignments(self,arglist):
        """Fetch and print alignments from nodes of the tree after treeMultipleAlign"""

        number=int(arglist[0])
        if len(arglist)>1:
            numofalign=int(arglist[1])
        else:
            numofalign=3

        self.alignment.getNodeAlignments(number,numofalign)
#        align=self.alignment.getNodeAlignments(number,numofalign)
#        self.show(Output.formatalign(align,self.seq))


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


            data=[x.strip().split() for x in open(filename).readlines() ]
            while(len(data[-1])==0):
                del data[-1]

            self.alignment=multiAlign.MultiAlignment(data,int(num_of_align),
                                                     float(Lambda), float(xi),
                                                     float(mu), float(nu),float(nuc_per_rotation))
            if not self.alignment:
                self.show("No multiple alignment for a reason or an other")
            else:
                self.show("Done")
                self.moreAlignments(int(num_of_align))
                self.show(str(len(self.alignment.bestAlignments)))
                self.show("Used time %g sec."%(self.alignment.secs_to_align))
                #for y in [(x.motif,x.score,zip(self.alignment.names,x.beginEnd,x.siteScore,x.siteSeqPos)) for x in self.alignment.bestAlignments[0]]:print y
                #print Output.formatalignGFF(self.alignment)
                #                print "goodAlign=",map(str,self.alignment.nextBest())
        except ValueError,e:
            raise CommandError("Error:",e)
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
            self.show(str(i))


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
            self.show("No multiple alignment to save!")
            return
        m="w"
        for i in self.malignment:
            i.setAlnLimit(minPairs)
            #if len(i)<2:continue
            if len([x for x in i.seqs if x in self.seq.getNames()])==len(i.seqs):
                i.strAln(self.seq)
            Output.savealign(str(i)+"\n",fname,m)
            m="a"

    def resetMatrices(self,arglist=None):
        "Arguments: none\nremoves all matrices"
        self.matlist=[]
        self.matdict={}
        self.resetTFBS()
        collect()

    def addMatrix(self, filenames):
        "adds new matrices from given files"
        self.resetTFBS()
        for f in filenames:
            try:
                m=Matrix.Matrix(f)
            except ValueError,e:
                raise CommandError("could not read matrix "+f)
            except IOError, (errno, strerror):
                raise CommandError("%s: %s" % (strerror, f))

            if self.matdict.has_key(f):
                pass
            else:
                self.show("adding %s, Info=%2.4g:"%(f,m.InfoContent))
                self.matlist.append(m)
                self.matdict[f]=m
        # Make the matrix names nicer.
        cpreflen=len(os.path.commonprefix([os.path.dirname(x.fname) for x in self.matlist]))
        try:
            cpreflen=len(os.path.dirname(self.matlist[0].fname[:cpreflen+1]))
            if cpreflen>0:
                cpreflen+=1
        except KeyError:
            pass

        for m in self.matlist:
            m.name=m.fname[cpreflen:]
        

            
    def printMatrices(self):
        "prints matrices to standard out"
        for i in range(len(self.matlist)):
            self.show("Matrix No %d (info=%g) %s"%(i,self.matlist[i].InfoContent,self.matlist[i].name))
            self.matlist[i].draw()

    def setPseudoCount(self,pseudoCnt):
        "Set the amount of pseudocounts on matricies."
        Matrix.setPseudoCount(pseudoCnt)
        self._initMatrixWeights()

    def setBGFreq(self,arglist=None):
        "Arguments: A C G T|bgSampleSequence\nBackground nucleotide frequencies. Removes markov background."
        if (not arglist==None):
            try:
                if len(arglist)==4:
                    tot=reduce(lambda x,y:float(x)+float(y),arglist,0.0)*1.0
                    self.A,self.C,self.G,self.T=map(lambda x:float(x)/tot,arglist)
                elif len(arglist)==1:
                    seq=self.seq.sequence(arglist[0]).getvalue()
                    self.A,self.C,self.G,self.T=[seq.count(x) for x in "ACGT"]
                    tot=float(self.A+self.C+self.G+self.T)
                    self.A,self.C,self.G,self.T=self.A/tot,self.C/tot,self.G/tot,self.T/tot
                else:
                    assert(1==0)
            except (ValueError,AssertionError,KeyError),e:
                raise CommandError("Invalid parameters as background frequencies.\nBackground distribution not set.",e)

        Matrix.setBGfreq(self.A,self.C,self.G,self.T)
        self._initMatrixWeights()

    def setMarkovBG(self,bgData,order=4):
        "Arguments: bgSampleSequence [order]\nBackground sample sequence and order of the model or saved background file."
        try:
            sampleStr=self.seq.sequence(bgData)
        except KeyError:
            filename=glob(bgData)
            sampleStr=eval(open(filename[0]).read())


        self.bg=_c_matrix.BackGround(sampleStr,int(order))

        Matrix.setMarkovBackground(self.bg)
        self._initMatrixWeights()


    def _initMatrixWeights(self):
        self.resetTFBS()
        map(Matrix.Matrix.initWeights,self.matlist)
            


    def printMatrixWeights(self):
        "prints matrix weights to standard out"
        for i in range(len(self.matlist)):
            try:
                bgStr="Order %d Markov"%(Matrix.Matrix.backGround.order)
            except AttributeError:
                bgStr="(%0.2f,%0.2f,%0.2f,%0.2f)"%(self.matlist[i].freqA,self.matlist[i].freqC,self.matlist[i].freqG,self.matlist[i].freqT)
            self.show("Matrix No %d %s bg=%s pCount=%g maxscore=%0.2f"%(i,self.matlist[i].name,bgStr,Matrix.Matrix.pseudoCount,self.matlist[i].maxscore))
            self.matlist[i].drawWeights()

    def removeMatrix(self, index):
        "removes matrix given by index"
        try:
            self.matlist.pop(index)
            self.resetTFBS()
        except IndexError:
            pass

    def resetSequences(self, arglist=None):
        "Arguments: none\nremoves all sequences"
        self.seq=Sequences()
        self.resetTFBS()
        collect()

    def addSequence(self, filenames):
        "adds sequences from files in FASTA format"
        for f in filenames:
            try:
                self.seq.addSequence(f)
            except Exception,e:
                raise CommandError("Can not read file %s"%(f),e)
            
    def removeSequence(self, name):
        "removes sequences given by name"
        removedSequences=self.seq.removeSequence(name)
        for n in removedSequences:
            try:
                del(self.__comp[n])
            except KeyError:
                pass
        collect()
    
    def getSeqNames(self):
        "returns the sequence names"
        return self.seq.getNames()


    def getTFBSAbsolute(self,cutoff=9.0):
        "computes the possible TFBS"
        Interface.getTFBS(self,cutoff,Matrix.CUTOFF_ABSOLUTE)

    def getTFBSpvalue(self,cutoff):
        "computes the possible TFBS"
        Interface.getTFBS(self,cutoff,Matrix.CUTOFF_PVALUE)
 
    def getTFBS(self, bound=0.1 , cutoffType=Matrix.CUTOFF_RELATIVE):
        "Scan for putative TFBS. For each matrix, score threshold is log2(bound) + maximum score."
        if hasattr(self,"tempFileName"):
            os.remove(self.tempFileName)
            del(self.tempFileName)

        #self.__comp={}
        if self.__comp_bound!=bound or self.__comp_cutoffType!=Matrix.CUTOFF_ABSOLUTE:
            self.resetTFBS()
        self.__comp_bound=bound
        self.__comp_cutoffType=cutoffType
        progress= 0.0
        totalMatches=0
        seqnumber=0
        matchSeqs=[name for name in self.seq.getNames() if not self.__comp.has_key(name)]
        for name in matchSeqs:
            seqnumber +=1
            self.show("Matching Sequence %s. %d of %d"%(name,seqnumber,len(matchSeqs)))
            self.__comp[name]={}
            try:
                startTime=time()
                self.__comp[name]=Matrix.getAllTFBS(self.seq.sequence(name),
                                                    bound,self.matlist,cutoffType)
                endTime=time()
            except (OverflowError,ValueError),e: # Zero, Negative
                raise CommandError("Need positive threshold!",e)
                return
#                self.__comp[name][m]=m.getTFBSbyRatio(self.seq.sequence(name),
#                                                      bound)
            try:
                foundMatches=sum([len(x) for x in self.__comp[name].values()])
                totalMatches+=foundMatches
                self.show("Found %d matches in %s\n"%(foundMatches,timeFormat(endTime-startTime)))
            except TypeError:
                self.show("name=",name)
                self.show("self.__comp=",self.__comp)
                #self.__gff=Output.get(self.__comp).split('\n')
            if totalMatches>100000:
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
            self.show("Storing temporary file",self.tempFileName)
        self.show("Writing temporary file")
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
            self.show("Storing file in gziped format.")
            if not filename[-2:]=="gz":
                filename=filename+".gz"
            try:
                shutil.copy(self.tempFileName,filename)
                os.remove(self.tempFileName)
                del(self.tempFileName)
            except IOError,e:
                self.show(e)
                filename=self.tempFileName
            return filename
        else:
            return Output.savematch(self.__comp, filename)

    def getmatchStr(self):
        "Get the list of matched TFBS:es as a string in GFF fromat"
        outStr=""

        if hasattr(self,"tempFileName"):
            try:
                tempFile=GzipFile(self.tempFileName,"r")
            except NameError:
                tempFile=open(self.tempFileName,"r")
            outStr=tempFile.read()
            tempFile.close()
        else:
            outStr=Output.get(self.__comp)
        return outStr

    def showmatch(self):
        "Prints the results to standard out"
        Output.showmatch(self.getmatchStr())


    def showMoreAlignments(self,count=1):
        """Gets more alignments but doesn't display them"""
        if not hasattr(self,"alignment"):
            return 
        if not self.alignment or not hasattr(self.alignment,"memSaveUsed"):
            self.show(self.alignment)
            return
        self.moreAlignments(count,self.alignment.nextBest)
        #print Output.formatalign(self.alignment,self.seq),

    def suboptimalAlignments(self,count=1):
        """Gets more alignments but doesn't display them"""
        try:
            self.moreAlignments(count,self.alignment.suboptimal)
        except (NotImplementedError,AttributeError):
            raise CommandError("Fetching suboptimal alignments is not supported.")
        #print Output.formatalign(self.alignment,self.seq),

    def suboptimalsDownTo(self,lowScore):
        """Gets alignments down to lowScore"""
        try:
            prevNumOfAlignments=-1
            while (len(self.alignment.bestAlignments)<1 or self.alignment.bestAlignments[-1][-1].score>lowScore) and len(self.alignment.bestAlignments)>prevNumOfAlignments:
                prevNumOfAlignments=len(self.alignment.bestAlignments)
                self.moreAlignments(1,self.alignment.suboptimal)
        except (NotImplementedError,AttributeError):
            self.show("Fetching suboptimal alignments is not supported.")

    def quit(self):
        "Exits the program"
        self.show("Exiting the program")
        sys.exit()


    def haveMatches(self):
        "Return true, if we have stored TFBS matches"
        return  len(self.__comp)>0 or hasattr(self,"tempFileName")
            

    def computeExpectationModel(self):
        """Compute the a and b for the E-value estimate E=exp(-(a+b*S))

        This is done by fitting a simple linear model
        log(r)=a+b*S+e
        where r is the rank of module, S is the score and e is an error
        term with expected value 0.

        Ten highest ranking modules are omited to account for 'real' signal"""

        #print "Called expectation model"
        if not (self.alignment and len(self.alignment.bestAlignments)>30):
            # Do nothing if there is no alignment data
            #print "Nothing to do",repr(self.alignment)
            #if self.alignment:
            #    print len(self.alignment.bestAlignments)
            return 
        Data=[(math.log(r+1),aS[-1].score) for r,aS in enumerate(self.alignment.bestAlignments)]

        # Do the regression on modules ranking 10-50
        Data=Data[10:50]

        n=1.0*len(Data)

        # Do the regression

        Sxy=sum([S*lr for (lr,S) in Data])
        Sx=sum([S for (lr,S) in Data])
        Sy=sum([lr for (lr,S) in Data])
        Sxx=sum([S**2 for (lr,S) in Data])
        try:
            self.alignment.beta=(n*Sxy-Sx*Sy)/(n*Sxx-Sx**2)
        except ZeroDivisionError:
            self.alignment.beta=1e99
        self.alignment.alpha=(Sy-self.alignment.beta*Sx)/n


        #self.beta=(n*sum([S*lr for (lr,S) in Data])-sum([S for (lr,S) in Data])*sum([lr for (lr,S) in Data]))/(n*sum([S**2 for (lr,S) in Data])-sum([S for (lr,S) in Data])**2)

        #self.alpha=(sum([lr for (lr,S) in Data])-self.alignment_A*sum([S for (lr,S) in Data]))/n

        # Do some diagnostics (R-squared i.e. the coefficient of determination)
        meanlr=Sy/n
        SSEnull=sum([(lr-meanlr)**2 for (lr,S) in Data])
        SSEmodel=sum([(lr-(self.alignment.alpha+self.alignment.beta*S))**2 for (lr,S) in Data])
        self.alignment.RMSE=math.sqrt(SSEmodel/n)
        try:
            self.alignment.Rsquared=1.0-SSEmodel/SSEnull
        except ZeroDivisionError:
            self.alignment.Rsquared=0.0
            
        

    def moreAlignments(self,num_of_align=1,fetcherFun=None):
        """Fetch more alignments from previously run alignment matrix"""

        if self.alignment.memSaveUsed==1:
            # Only nextBest fetcher is possible if memSave was used.
            fetcherFun=self.alignment.nextBest
            
        for i in range(num_of_align):
            if self.alignment.memSaveUsed==1 and self.alignment.askedResults<=len(self.alignment.bestAlignments):
                self.show("Can't give more alignments. Don't remember those")
                break
            if not fetcherFun:
                goodAlign=self.alignment.suboptimal()
            else:
                goodAlign=fetcherFun()
            if not goodAlign:
                break
            goodAlign.reverse() # For old times sake
            #               0 1 2     3      4[0]   4[1]   5[0]   5[1]  6
            # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
            self.alignment.bestAlignments.append(goodAlign)
        self.computeExpectationModel()

    def align(self, filename='.', num_of_align=3,
              Lambda=1.0, xi=1.0, mu=0.5,nu=1.0, nuc_per_rotation=10.4,
              firstSeq=None,secondSeq=None):
        "aligns the computed BS or some given in file"
        if hasattr(self,"alignment"):
            del(self.alignment)

        if filename!='.':
            files=glob(filename)
            if len(files)==0:
                raise AttributeError("Can't find input file!")
            elif len(files)>1:
                raise AttributeError("Ambiguous input file glob. Could be "+",".join(files))
            filename=files[0]
            
        if filename=='.' and  hasattr(self,"tempFileName"):
            filename=self.tempFileName

        collect()

        if filename=='.':
            if len(self.__comp)==0:
                return "No binding sites"
            else:
                self.alignment=align.aligndata(Output.get(self.__comp),num_of_align,
                                               Lambda, xi,
                                               mu, nu,nuc_per_rotation,firstSeq,secondSeq)

        else:
            self.alignment= align.alignfile(filename, num_of_align,Lambda,
                                            xi, mu, nu,nuc_per_rotation,firstSeq,secondSeq)
        self.moreAlignments(num_of_align)
        if self.alignment:
            #Interface.showalignSTDO(self)
            assert self.alignment.x_name!=self.alignment.y_name
            self.show("Used time %g sec."%(self.alignment.secs_to_align))
            return 1
        else:
            return None

    def randomize_full(self,arglist):
        """Arguments: none\nFully randomize the list of matched TFBS:es."""
        if len(self.__comp)==0:
            return

        import random
        # self.__comp = {'Sequence':{Matrix:{(position,strand,[ambig,allele,snpPos,scoreDif]):(Score,altScore)}}}

        self.__rcomp={}
        for sequence,matrixDict in self.__comp.items():
            pos=[]
            for matrix,hitDict in matrixDict.items():
                
                pos.extend([position for position,strand,SNPlist in hitDict.keys()])
            random.shuffle(pos)
            for matrix,hitDict in matrixDict.items():
                randHitDict=dict([((pos.pop(),strand,SNPlist),score) for ((position,strand,SNPlist),score) in hitDict.items()])
                matrixDict[matrix]=randHitDict
        
        self.show("Fully, Independently Randomized. Preserving positions.")


    def savealignGFF(self,filename=""):
        "Saves the results in GFF format"
        if 1:
#	try:
            return Output.savealign(Output.formatMultiAlignGFF(self.alignment,self.seq), filename)
#        except AttributeError:
#            raise CommandError("No alignment to save")

    def savealignAnchor(self,filename=""):
        "Saves the results in Anchor format"
        try:
            return Output.savealign(Output.formatalignCHAOS(self.alignment), filename)
        except AttributeError:
            raise CommandError("No alignment to save")


    def savealign(self, filename=''):
        "Saves the results"
        try:
            return Output.savealign(Output.formatalign(self.alignment,self.seq), filename)
        except AttributeError:
            raise CommandError("No alignment to save")


    def showalignSTDO(self):
        "Prints the alignment to standart out"
        if hasattr(self,"alignment"):
            self.show(Output.formatalign(self.alignment,self.seq))
    
    def savepwbase(self, filename=''):
        "Saves the pairwise alignments the multiple alignment was based on"
        try:
            return Output.savealign(Output.formatpwbase(self.alignment), filename)
        except AttributeError:
            raise CommandError("No alignment to save")


    def showpwbaseSTDO(self):
        "Prints the pairwise alignments the multiple alignment was based on to standart out"
        if hasattr(self,"alignment"):
            self.show(Output.formatpwbase(self.alignment))

    def about(self):
        "Information about the program"
        return """Enhancer Element Locator EEL
Version $Version$
Contact: Kimmo Palin (kimmo.palin@helsinki.fi)

Big thanks to Matthias Berg (bergm@web.de) for initial coding.
"""
#         Uni 18, Zi. 2220
#         D-66123 Saarbruecken
#         Germany


