"""Support for commandline interface.

This module takes care of the commandline history and parsing.
"""

from time import localtime
from glob import glob
from Interface import Interface
from popen2 import popen3
import re
import string


# $Log$
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
    histfile = os.path.join(os.environ["HOME"], ".mabshist")
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




class Commandline(Interface):
    """This is a command line based user interface.

    The mabs commands and most of the default values are set here."""
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
                         "setMarkovBG":        (self.setMarkovBG,1),
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
                         'more':               (self.moreAlignment,0)}

    def run(self):
        "waits for std input and executes these commands"
        print 'Type "help" for more information'
        while(1):
            try:
                token=[]
                while len(token)==0:
                    #print "> ",
                    try:
                        # read tokens for std input
                        token=raw_input("> ").split()
                    #except EOFError:
                    #    self.quit([])
                    except KeyboardInterrupt:
                        print

                # if command exists...
                if self.isCommand(token[0]): 
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
                        except StandardError,e:
                            print "#"*67
                            print "Software error encountered! Please email this"
                            print "error message to kimmo.palin@helsinki.fi"
                            print "Error while processing command:",token
                            import traceback,sys
                            traceback.print_exc(file=sys.stdout)
                            print "#"*67
                else:
                    print token[0],": command not found"
            except KeyboardInterrupt:
                print "Keyboard Interrupt: Aborted"

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
                    except StandardError,e:
                        print "#"*67
                        print "Software error encountered! Please email this"
                        print "error message to kimmo.palin@helsinki.fi"
                        print "Error while processing command:",token
                        import traceback,sys
                        traceback.print_exc(file=sys.stdout)
                        print "#"*67
            else:
                print token[0],": command not found"
        except KeyboardInterrupt:
            print "Keyboard Interrupt: Aborted"
        except IndexError:
            pass

    def isCommand(self, command):
        return self.__commands.has_key(command)


    def setBGFreq(self,arglist=None):
        "Arguments: A C G T\nBackground nucleotide frequencies. Removes markov background."
        if not arglist==None:
            tot=reduce(lambda x,y:float(x)+float(y),arglist,0.0)*1.0
            self.A,self.C,self.G,self.T=map(lambda x:float(x)/tot,arglist)
            
        for m in self.matlist:
            m.setBGfreq(self.A,self.C,self.G,self.T)

    def saveMarkovBackground(self,arglist=None):
        "Arguments: filename\nName of the file where to store the background model.\n"
        grams=None
        try:
            grams=self.bg.giveGramVector()
        except AttributeError,e:
            print e
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

            

    def setMarkovBG(self,arglist=None):
        "Arguments: bgSampleSequence [order]\nBackground sample sequence and order of the model or saved background file."
        try:
            sampleStr=self.seq.sequence(arglist[0])
        except KeyError:
            filename=glob(arglist[0])
            try:
                sampleStr=eval(open(filename[0]).read())
            except Exception,e:
                print "No such sequence or saved background (%s). Markov background not set."%(str(e))
                return
                pass
        order=4
        if len(arglist)>1:
            order=int(arglist[1])

        import matrix
        self.bg=matrix.BackGround(sampleStr,order)

        for m in self.matlist:
            m.setMarkovBackground(self.bg)
            

    def setPseudoCnt(self,arglist=("1.0")):
        "Arguments: [pseudocount]\nSet the amount of pseudocounts on matricies. Default 1.0"
        for m in self.matlist:
            m.setPseudoCount(string.atof(arglist[0]))

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
                self.setBGFreq()
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
        """Arguments: filelist
        Gzipped and Fasta formated sequence files. One sequence in file."""
        for filestring in arglist:
            filenames=glob(filestring)
            for filename in filenames:
                self.seq.addSingleSequence(filename)
            if not filenames:
                print "file not found:", filestring
        
    def addSequence(self, arglist):
        "Arguments: filelist\nreads sequences from files"
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
                Interface.addSequence(self, filenames)
            else:
                print "file not found:", filestring

    def printSeqNames(self, arglist):
        "Arguments: none\nprints the names of the sequences"
        print str(self.seq),
        #names= Interface.getSeqNames(self)
        #names.sort()
        #for n in names:
        #    print n

    def removeSequence(self, arglist):
        "Arguments: Sequencename\nremoves a sequence"
        seqname=arglist[0]
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
            print "arglist:",repr(arglist)
            print arglist[0],"is not a number"


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
            print "arglist:",repr(arglist)
            print arglist[0],"is not a number"
            
    def showmatch(self,arglist):
        "Arguments: none\nprints the computed scores to stdout"
        Interface.showmatch(self)

    def help(self, arglist):
        "Arguments: none\nprints this help"
        minus=''
        if arglist==['add -']: minus='-'
        values=self.__commands.values()
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
The default filename is 'mabs_[Date+Time].gff'
e.g. mabs_2003_8_27_15_48.gff"""
        filename=''
        if len(arglist):
            filename=arglist[0]
        filename = Interface.savematch(self,filename)
        if filename:
            print"results saved in", filename

    def resetMatrices(self, arglist):
        "Arguments: none\nremoves all matrices"
        Interface.resetMatrices(self)

    def resetSequences(self, arglist):
        "Arguments: none\nremoves all sequences"
        Interface.resetSequences(self)

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

            if not Interface.align(self, filename, int(num_of_align),
                            float(Lambda), float(xi),
                            float(mu), float(nu),float(nuc_per_rotation)):
                print "No alignment for a reason or an other"
            
        except ValueError:
            print "Error: unallowed arguments passed to 'align'"


    def moreAlignment(self,arglist=(1,)):
        """Arguments: [number of alignments]
        prints more alignments from previously run alignment matrix"""
        if len(arglist)>0:
            count=string.atoi(arglist[0])
        else:
            count=1
        self.showMoreAlignments(count)



    def getBaseSaveName(self):
        """Assistant function that returns the default basename for output files"""
        a=localtime()
        filename='mabs_'+str(a.tm_year)+'_'+str(a.tm_mon)+'_'+str(a.tm_mday)+'_'+str(a.tm_hour)+'_'+str(a.tm_min)
        return filename

        
    def savealign(self, arglist):
        """Arguments: [filename]
saves the alignment to disk
The default filename is 'mabs_[Date+Time].align'
e.g. mabs_2003_9_16_11_48.align"""
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
The default filename is 'mabs_[Date+Time]_align.gff'
e.g. mabs_2003_9_16_11_48_align.gff"""
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
The default filename is 'mabs_[Date+Time]_align.anc'
e.g. mabs_2003_9_16_11_48_align.anc"""
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
