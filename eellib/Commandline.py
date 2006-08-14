# -*- coding: UTF-8 -*-

"""Support for commandline interface.

This module takes care of the commandline history and parsing.
"""

from time import localtime
from glob import glob
from Interface import Interface
from popen2 import popen3
import re
import string

import _c_matrix


# $Log$
# Revision 1.21  2006/04/05 08:30:07  kpalin
# Regular expressions for sequence removal and commands for multiple alignment.
#
# Revision 1.20  2005/05/19 07:49:35  kpalin
# Merged Waterman-Eggert style suboptimal alignments and
# SNP matching.
#
# Revision 1.19.2.3  2005/05/09 07:14:57  kpalin
# Matrix distances with hidden command.
#
# Revision 1.19.2.2  2005/04/12 09:09:57  kpalin
# Now help accepts parameters and ps outputs sequence description along
# with names.
#
# Revision 1.19.2.1  2005/03/31 13:30:54  kpalin
# Added command 'suboptimal' which is like 'more' but gives
# real suboptimal results instead of next best from the alignment matrix.
#
# Revision 1.19  2005/03/22 13:17:13  kpalin
# Merged some fixes surfacing from testing the public version.
#
# Revision 1.18  2005/03/09 07:40:03  kpalin
# Fixed order setting in setMarkovBG
#
# Revision 1.17  2005/03/08 10:39:16  kpalin
# Fixed Matrix add/Markov background setting.
#
# Revision 1.16  2005/01/13 13:16:42  kpalin
# Moved the requesting of sequences to be aligned to Python side
# of stuff. Much better.
#
# Revision 1.15  2005/01/13 09:38:30  kpalin
# Added catch for EOFError in no-gui
#
# Revision 1.14  2005/01/12 13:34:55  kpalin
# Added Tkinter/Tix Graphical user interface and command -no-gui to
# avoid it.
#
# Revision 1.13  2005/01/07 13:41:25  kpalin
# Works with py2exe. (windows executables)
#
# Revision 1.12  2004/12/22 11:14:24  kpalin
# Some fixes for better distributability
#
# Revision 1.11  2004/12/22 08:02:59  kpalin
# Hopefully more IO efficient TFBS search.
#
# Revision 1.10  2004/12/14 13:08:05  kpalin
#
# Name change from MABS to EEL (Enhancer Element Locator / Monty Python pun
# "My hovercraft is full of EELs" )
#
# Revision 1.9  2004/07/30 12:21:08  kpalin
# Commands for multiple alignment and discouragement of using greedy
# alignment.
#
# Revision 1.8  2004/04/08 13:05:38  kpalin
# Some cleaning and fixes.
#
# Revision 1.7  2004/03/03 09:27:25  kpalin
# Added interface for multiple alignment.
#
# Revision 1.6  2004/02/26 11:59:49  kpalin
# Added gracefull error handling for mabs. Now it doesn't exit after
# encountering an error but report it and tells to send the error
# message to me.
#
# Revision 1.5  2004/01/28 08:48:32  kpalin
# Updated docstrings
#
# Revision 1.4  2004/01/14 10:05:57  kpalin
# Generated documentation
#
# Revision 1.3  2004/01/13 07:55:01  kpalin
# Paljon kaikenlaista.
#
# Mahdotonta muistaa kaikkea.
#
# Revision 1.2  2003/12/29 12:43:17  kpalin
# Interface class repaired to enable alignment from gzip:ed temporary files.
#
# Ilmeisesti jotain uutta. En tiedä mitä.
#


try:
    # Use historyfile
    import os
    import readline
    histfile = os.path.join(os.environ["HOME"], ".eelhist")
    try:
        readline.read_history_file(histfile)
    except IOError:
        pass
    import atexit
    atexit.register(readline.write_history_file, histfile)
    del os, histfile
except ImportError:
    print "No command history available."
except Exception:
    pass

try:
    import os
except ImportError:
    pass





class Commandline(Interface):
    """This is a command line based user interface.

    The EEL commands and most of the default values are set here."""
    def __init__(self):
        Interface.__init__(self)
        self.A,self.C,self.G,self.T=0.25,0.25,0.25,0.25
        
        # This is a map from commands to their functions and
        # at least required number of arguments
        self.__commands={'about':              (self.about,0),
                         'h':                  (self.help,0),
                         'help':               (self.help,0),
                         '?':                  (self.help,0),
                         'addMatrix':          (self.addMatrix,1),
                         'am':                 (self.addMatrix,1),
                         'removeMatrix':       (self.removeMatrix,1),
                         'rm':                 (self.removeMatrix,1),
                         'align':              (self.align,0),
                         'resetMatrices':      (self.resetMatrices,0),
                         'resm':               (self.resetMatrices,0),
                         'printMatrices':      (self.printMatrices,0),
                         'pm':                 (self.printMatrices,0),
                         'pmw':                (self.printMatrixWeights,0),
                         'printMatrixWeights': (self.printMatrixWeights,0),
                         'addSequence':        (self.addSequence,1),
                         'as':                 (self.addSequence,1), 
                         'removeSequence':     (self.removeSequence,1),
                         'rs':                 (self.removeSequence,1),
                         'resetSequences':     (self.resetSequences,0),
                         'ress':               (self.resetSequences,0),
                         'reset':              (self.reset,0),
                         "setBGfreq":          (self.setBGFreq,4),
                         "setMarkovBG":        (self._setMarkovBG,1),
                         'printSeqNames':      (self.printSeqNames,0),
                         'ps':                 (self.printSeqNames,0),
                         'getTFBS':            (self.getTFBS,0),
                         'getTFBSabsolute':    (self.getTFBSabsolute,0),
                         'showmatch':          (self.showmatch,0),
                         'sm':                 (self.showmatch,0),
                         'q':                  (self.quit,0),
                         'quit':               (self.quit,0),
                         'savematch':          (self.savematch,0),
                         'savealign':          (self.savealign,0),
                         'savealignGFF':       (self.savealignGFF,0),
                         'savealignAnchor':       (self.savealignAnchor,0),
                         'showalign':          (self.showalign,0),
                         'sa':                 (self.showalign,0),
                         'setpseudocount':     (self.setPseudoCnt,1),
                         'addSingleSequence':      (self.addSingleSequence,1),
                         'ass':               (self.addSingleSequence,1),
                         'saveMarkovBackground':  (self.saveMarkovBackground,1),
                         'more':               (self.moreAlignment,0),
                         'suboptimals':               (self.suboptimalAlignment,0),
                         '__multipleAlignGreedy':      (self.multiAlignGreedy,1),
                         'multipleAlign':      (self.multiAlign,0),
                         '__showMultiAlign':     (self.showMultiAlign,0),
                         '__saveMultiAlign':     (self.saveMultiAlign,1),
                         'no-gui':               (self.no_gui,0),
                         '__computeKLdistances': (self.showKLdist,1),
                         '__computeEscores': (self.showExpectedScores,1),
                         '__head': (self.getHead,1),
                         '__randomize_full': (self.randomize_full,0)
                         }
        # Add directory commands if available.
        if globals().has_key("os") and hasattr(os,"getcwd") and hasattr(os,"chdir"):
            self.__commands.update({'dir':                  (self.dirlist,0),
                                    'cd':                   (self.chgdir,0)})

    def run(self):
        "waits for std input and executes these commands"
        print 'Type "help" for more information'
        while(1):
            token=[]
            while len(token)==0:
                try:
                    # read tokens for std input
                    token=raw_input("> ")
                #except EOFError:
                #    self.quit([])
                except KeyboardInterrupt:
                    print


            self.execute(token)

    def execute(self, command):
        "executes the given command (like 'run' does)"
        try:
            token=command.split()
            
            # if command exists...
            if len(token) and self.__commands.has_key(token[0]):
                comm=self.__commands[token[0]][0]
                num_of_args=self.__commands[token[0]][1]

                # if command gets too less arguments...
                if num_of_args > len(token)-1:
                    s=''
                    if num_of_args!=1: s='s'
                    print token[0], "needs at least ",\
                          num_of_args, "argument"+s
                # else execute command with its argumments
                else:
                    try:
                        comm(token[1:])
                    except KeyboardInterrupt:
                        raise
                    except StandardError,e:
                        import sys
                        print "#"*67
                        print "Software error encountered! Please email this"
                        print "error message to kimmo.palin@helsinki.fi"
                        print "Platform:",sys.platform
                        print "Version: '%s'"%(sys.version)
                        print "Error while processing command:",token
                        import traceback,sys
                        traceback.print_exc(file=sys.stdout)
                        try:
                            print "API-version: %s"%(str(sys.api_version))
                        except AttributeError:
                            pass
                        print "#"*67
            else:
                print token[0],": command not found"
        except KeyboardInterrupt:
            print "Keyboard Interrupt: Aborted"
        except IndexError:
            pass

    def isCommand(self, command):
        return self.__commands.has_key(command)



    def getHead(self,arglist):
        "Arguments: sites\nCut site sequences to length sites"
        self.head(int(arglist[0]))
        
    def no_gui(self,arglist=None):
        "Arguments: none\nGives command line interface"
        try:
            self.run()
        except EOFError:
            pass



    def dirlist(self,arglist=None):
        "Arguments: [path/pattern]\nList files matching the given pattern. Defaults to \"*\" in current working directory.\n"
        try:
            dirPath=arglist[0]
        except IndexError:
            dirPath="*"

        
        files=glob(dirPath)
        self.showFileList(files)


    def chgdir(self,arglist=None):
        "Arguments: [path]\nChange or display the current working directory."
        try:
            if arglist:
                os.chdir(arglist[0])
        
            newPath=os.getcwd()
            self.show(newPath)
        except OSError,e:
            self.show(str(e))
            
            


    def saveMarkovBackground(self,arglist=None):
        "Arguments: filename\nName of the file where to store the background model.\n"
        grams=None
        try:
            grams=self.bg.giveGramVector()
        except AttributeError,e:
            print "No markov background to save"
            return

        if grams:
            try:
                fout=open(arglist[0],"w")
                fout.write(repr(grams))
                fout.close()
            except IOError,e:
                print "Could not save"
                print e
                return

            

    def _setMarkovBG(self,arglist=None):
        "Arguments: bgSampleSequence [order]\nBackground sample sequence and order of the model or saved background file."
        try:
            order=int(arglist[1])
        except IndexError:
            order=4
        except ValueError:
            print "Invalid order. Markov background not set."
            return
        
        try:
            datName=arglist[0]
        except IndexError:
            print "Missing file or sample sequence name. Markov background not set."
            return
        try:
            self.setMarkovBG(arglist[0],order)
        except Exception,e:
            print "No such sequence or saved background (%s). Markov background not set."%(str(e))
            return


    def setPseudoCnt(self,arglist=("1.0")):
        "Arguments: [pseudocount]\nSet the amount of pseudocounts on matricies. Default 1.0"
        try:
            Interface.setPseudoCount(self,float(arglist[0]))
        except ValueError:
            print "Invalid pseudocount. Pseudocount not set."
            return

    # The following functions are executable from the command line (see 'run')
    # The Argument 'arglist' is a list of the command line arguments
    # The '__doc__'-strings are used by the 'help' function, so pay attention
    # when editing these.

    def addMatrix(self, arglist):
        "Arguments: filelist\nreads matrices from files"
        for filestring in arglist:
##            if (filestring.find("/")==-1):
##                filestring = "./" + filestring
##            # dividing the fiestring in path and filename
##            filematch= re.match("(.*)/(.*)", filestring)
##            # using the 'find' command to allow things like '*' in filename
##            syscomm=('find '+filematch.group(1)+' -maxdepth 1 -xtype f -name "'
##                     +filematch.group(2)+'"')
##            filenames= popen3(syscomm)[0].read().split()
            filenames=glob(filestring)
            if filenames:
                Interface.addMatrix(self, filenames)
            else:
                print "file not found:", filestring

    def removeMatrix(self, arglist):
        "Arguments: Matrixnumber\nremoves a matrix"
        try:
            index=string.atoi(arglist[0])
            Interface.removeMatrix(self, index)
        except ValueError:
            print arglist[0],"is not a number"

    def printMatrixWeights(self,arglist):
        "Arguments: none\nprints the matrix weights (with background)"
        Interface.printMatrixWeights(self)
        
    def printMatrices(self, arglist):
        "Arguments: none\nprints the matrices"
        Interface.printMatrices(self)


    def addSingleSequence(self,arglist):
        "Arguments: filelist\nGzipped and Fasta formated sequence files. One sequence in file."
        for filestring in arglist:
            filenames=glob(filestring)
            for filename in filenames:
                self.seq.addSingleSequence(filename)
            if not filenames:
                print "file not found:", filestring
        
    def addSequence(self, arglist):
        "Arguments: filelist\nreads sequences from files"
        for filestring in arglist:
            filenames=glob(filestring)
            if filenames:
                Interface.addSequence(self, filenames)
            else:
                print "file not found:", filestring

    def printSeqNames(self, arglist):
        "Arguments: none\nprints the names of the sequences"
        print str(self.seq)
        #names= Interface.getSeqNames(self)
        #names.sort()
        #for n in names:
        #    print n

    def removeSequence(self, arglist):
        "Arguments: Sequencename\nremoves a sequence"
        for seqname in arglist:
            Interface.removeSequence(self, seqname)
        

    def getTFBS(self, arglist):
        """Arguments: [bound]
computes the scores of all matrices and all sequences which are
better than bound*maxscore. maxscore is the highest reachable
score of the actual matrix with respect to the background
The default value for bound is 0.1"""
        bound=0.1
        try:
            if(len(arglist)>0):
                bound= string.atof(arglist[0])
            Interface.getTFBS(self, bound)
        except ValueError:
            print "getTFBS requires an numeric argument!"


    def getTFBSabsolute(self, arglist):
        """Arguments: [cutoff]
computes the scores of all matrices and all sequences which are
better than cutoff.
The default value for cutoff is 9.0"""
        cutoff=9.0
        try:
            if(len(arglist)>0):
                cutoff= string.atof(arglist[0])
            Interface.getTFBSAbsolute(self, cutoff)
        except ValueError:
            print "getTFBSabsolute requires an numeric argument!"
            


    def showKLdist(self,arglist):
        "Arguments: File to save the result\nComputes the distance matrix between all pairs of TFBS matrixes."
        KL=dict([(m.name,{}) for m in self.matlist])
        try:
            fout=open(arglist[0],"w")
        except IOEerror,e:
            print "Could not save %s:"%(arglist[0]),e
            return

        for m in self.matlist:
            for n in self.matlist:
                KL[m.name][n.name]=m.minimumKLdistance(n)

        names=KL.keys()
        names.sort()
        fout.write("\t".join(["names"]+names)+"\n")
        fout.writelines(["%s\t%s\n"%(name,"\t".join([str(KL[name][n][0]) for n in names])) for name in names])


    def showExpectedScores(self,arglist):
        "Arguments: File to save the result\nComputes the expected scores of matrices, given the other."
        KL=dict([(m.name,{}) for m in self.matlist])
        try:
            fout=open(arglist[0],"w")
        except IOEerror,e:
            print "Could not save %s:"%(arglist[0]),e
            return

        for m in self.matlist:
            for n in self.matlist:
                KL[m.name][n.name]=m.maximumExpectedScore(n)

        names=KL.keys()
        names.sort()
        fout.write("\t".join(["RowGivenCol"]+names)+"\n")
        fout.writelines(["%s\t%s\n"%(name,"\t".join([str(KL[name][n][0]) for n in names])) for name in names])
        
    def showmatch(self,arglist):
        "Arguments: none\nprints the computed scores to stdout"
        Interface.showmatch(self)

    def help(self, arglist):
        "Arguments: none\nprints this help"
        minus=''
        if arglist==['add -']:
            minus='-'
            arglist==[]
        if "full" in arglist:
            values=self.__commands.values()
        else:
            values=[x[1] for x in self.__commands.items() if x[0][:2]!="__" and (len(arglist)==0 or x[0] in arglist)]
        values.sort()
        for v in values:
            if values.count(v)>1:
                values.remove(v)

        helplist=[]
        for v in values:
            keys=[]
            for k in self.__commands.keys():
                if self.__commands[k]==v:
                    keys.append(minus+k)
            keys.sort()
            entry= str(keys)[1:-1]
            for line in v[0].__doc__.split('\n'):
                entry+= '\n\t'+line
            helplist.append(entry)

        helplist.sort()
        for entry in helplist:
            print entry
            print
            
                
    def quit(self, arglist):
        "Arguments: none\nto exit the program"
        Interface.quit(self)

    def savematch(self, arglist):
        """Arguments: [filename]
saves the results of the matching in gff format
See http://www.sanger.ac.uk/Software/formats/GFF/
The default filename is 'eel_[Date+Time].gff'
e.g. eel_2003_8_27_15_48.gff"""
        filename=''
        if len(arglist):
            filename=arglist[0]
        filename = Interface.savematch(self,filename)
        if filename:
            print"results saved in", filename



    def reset(self,arglist):
        "Arguments: none\nremoves all matrices and sequences"
        Interface.resetMatrices(self)
        Interface.resetSequences(self)

        
    def align(self, arglist):
        """Arguments: [filename[,num_of_align,[lambda[,xi[,mu[,nu,[,nuc_per_rotation]]]]]]]
aligns the computed BS or optional the BS from a gff file
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

            try:
                if not Interface.align(self, filename, int(num_of_align),
                                       float(Lambda), float(xi),
                                       float(mu), float(nu),float(nuc_per_rotation)):
                    print "No alignment for a reason or an other"
            except AttributeError,val:
                if len(val.args)==1 or (len(val.args)>1 and len(val.args[1])<2):
                    print val.args[0]
                else:
                    msg,seqs=val
                    # This function is expected to be defined in subclass
                    firstSeq,secondSeq=self.selectTwoSequences(msg,seqs)
                    if not Interface.align(self, filename, int(num_of_align),
                                           float(Lambda), float(xi),
                                           float(mu), float(nu),
                                           float(nuc_per_rotation),
                                           firstSeq,secondSeq):
                        print "No alignment for a reson or an another"
        except ValueError,e:
            print "Error: unallowed arguments passed to 'align'",e


    def selectTwoSequences(self,message,sequences):
        print message
        while 1:
            print "Select first sequence"
            seqs=zip(range(len(sequences)),sequences)

            for i,seq in seqs:
                print "%d %s"%(i,seq)

            print
            firstSeqCode=raw_input("Give Sequence: ")
            try:
                firstSeqCode=int(firstSeqCode)
                assert(0<=firstSeqCode<len(sequences))
                firstSeq=sequences[firstSeqCode]
                del(sequences[firstSeqCode])
                break
            except(AssertionError,ValueError):
                pass

        while 1:
            print "Select second sequence"
            seqs=zip(range(len(sequences)),sequences)

            for i,seq in seqs:
                print "%d %s"%(i,seq)

            print
            secondSeqCode=raw_input("Give Sequence: ")
            try:
                secondSeqCode=int(secondSeqCode)
                assert(0<=secondSeqCode<len(sequences))
                secondSeq=sequences[secondSeqCode]
                del(sequences[secondSeqCode])
                break
            except(AssertionError,ValueError):
                pass

        return (firstSeq,secondSeq)


    def moreAlignment(self,arglist=(1,)):
        """Arguments: [number of alignments]
        prints more alignments from previously run alignment matrix"""
        if len(arglist)>0:
            count=string.atoi(arglist[0])
        else:
            count=1
        self.showMoreAlignments(count)

    def suboptimalAlignment(self,arglist):
        """Arguments: [number of alignments]
        prints more suboptimal alignments from previously run matrix"""
        
        if len(arglist)>0:
            count=string.atoi(arglist[0])
        else:
            count=1
        self.suboptimalAlignments(count)


    def getBaseSaveName(self):
        """Assistant function that returns the default basename for output files"""
        a=localtime()
        filename='eel_'+str(a.tm_year)+'_'+str(a.tm_mon)+'_'+str(a.tm_mday)+'_'+str(a.tm_hour)+'_'+str(a.tm_min)
        return filename

        
    def savealign(self, arglist):
        """Arguments: [filename]
saves the alignment to disk
The default filename is 'eel_[Date+Time].align'
e.g. eel_2003_9_16_11_48.align"""
        filename=''
        if len(arglist):
            filename=arglist[0]
        else:
            filename=self.getBaseSaveName()+'.align'
        filename = Interface.savealign(self,filename)
        if filename:
            print"results saved in", filename


    def savealignGFF(self, arglist):
        """Arguments: [filename]
saves the alignment to disk in GFF format
The default filename is 'eel_[Date+Time]_align.gff'
e.g. eel_2003_9_16_11_48_align.gff"""
        filename=''
        if len(arglist):
            filename=arglist[0]
        else:
            filename=self.getBaseSaveName()+'.align.gff'
            
        filename = Interface.savealignGFF(self,filename)
        if filename:
            print"results saved in", filename


    def savealignAnchor(self, arglist):
        """Arguments: [filename]
saves the alignment to disk in Anchor format for DIALIGN
The default filename is 'eel_[Date+Time]_align.anc'
e.g. eel_2003_9_16_11_48_align.anc"""
        filename=''
        if len(arglist):
            filename=arglist[0]
        else:
            filename=self.getBaseSaveName()+'.align.anc'
            
        filename = Interface.savealignAnchor(self,filename)
        if filename:
            print"results saved in", filename

    def showalign(self, arglist):
        "Arguments: none\nprints the computed alignment to stdout"
        self.showalignSTDO()    
        
    def about(self, arglist):
        "Arguments: none\nPrints Information about the program"
        print Interface.about(self)
