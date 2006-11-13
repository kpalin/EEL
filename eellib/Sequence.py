"""Sequence API.

Fasta file input is handled here."""
import string
from cStringIO import StringIO
#from xreadlines import xreadlines

import sys
if sys.platform!='win32':
    try:
        from gzip import GzipFile
    except ImportError:
        print "No gzip available."


#
# $Log$
# Revision 1.11  2006/04/05 08:30:07  kpalin
# Regular expressions for sequence removal and commands for multiple alignment.
#
# Revision 1.10  2005/07/08 07:56:43  kpalin
# Fixed a problem with seek() and multiple open()s in addSequence().
# Now it should work better with pipes.
#
# Revision 1.9  2005/05/19 07:49:35  kpalin
# Merged Waterman-Eggert style suboptimal alignments and
# SNP matching.
#
# Revision 1.8.2.1  2005/04/12 09:11:19  kpalin
# Handles also the descriptions along with sequence names.
#
# Revision 1.8  2005/03/03 09:07:05  kpalin
# Added handling for EMBL formated sequences. Not just FASTA.
#
# Revision 1.7  2005/01/07 13:41:25  kpalin
# Works with py2exe. (windows executables)
#
# Revision 1.6  2004/12/17 12:17:21  kpalin
# Fixed depreciated xreadline
#
# Revision 1.5  2004/03/03 09:15:43  kpalin
# Corrected a bug from SingleSequence.getslice() and added better
# resembelance to mapping object.
#
#
    


class SingleSequence:
    """Interface for single sequence files.

    This is useful for very large sequences and interface works more or less like a string."""
    def __init__(self,file):
        """Initialization.

        Parameter file is an open file like object"""
        self.file=file

        self.start=self.file.tell()
        self.file.seek(0,2)

        self.end=self.file.tell()
        self.file.seek(self.start,0)

        self.len=self.end-self.start
        
    def __getattr__(self,attr):
        """Attribute override.

        len() mehod is overriden"""
        if attr=="__len__":
            return self.__len__
        
        return getattr(self.file,attr)

    def __hasattr__(self,attr):
        return getattr(self.file,attr)

    def __getitem__(self,i):
        return self.__getslice__(i,i+1)

    def __getslice__(self,i,j):
        #print "Getting slice [%d:%d]"%(i,j)
        self.file.seek(i+self.start)
        c=self.file.read(j-i)
        self.file.seek(self.start)

        return c
    
    def __len__(self):
        return self.len


class UnsupportedTypeException(Exception):
    def __init__(self,v="Unsupported sequence format"):
        Exception.__init__(self,v)

class Sequences:
    "represents DNA-Sequences"
    def __init__(self,filename=0):
        "reads seqences from file"
        self.__Seq={}
        self.__Desc={}
        if filename:
            self.addSequence(filename)



    def resetSeqPositions(self):
        "Helper function. Resets the sequence position after every read of the sequence"
        for seq,pos in self.__Seq.values():
            seq.seek(pos)

    def addSingleSequence(self,filename):
        """Adds a single sequence from fasta formated and gziped file.
        Useful for huge sequences, e.g. whole chromosomes"""
        try:
            try:
                File=GzipFile(filename,'r')
                line=File.readline()
            except (NameError,IOError):
                File=open(filename,"r")
                line=File.readline()
            

            #print line
            if line[0]==">":
                name=string.split(line[1:].strip())[0]
                self.__Seq[name]=(SingleSequence(File),File.tell())
                print filename, "added"
            else:
                raise KeyError
        except IOError, e:
            if type(e)==type((1,2)):
                (errno, strerror)=e
                print "I/O error(%s): %s" % (errno, strerror)
            else:
                print e
        except KeyError:
            print filename,"in not in FASTA format!"



    def readFasta(self,fastaStr):
        toGetValue={}
        for line in fastaStr.split("\n"):
            line=line.strip()
            if line:
                if line[0]=='>':
                    parts=line[1:].strip().split(" ",1)
                    try:
                        name,description=parts
                    except ValueError:
                        description=""
                        name=parts[0]
                    self.__Desc[name]=description
                    toGetValue[name]=1
                    self.__Seq[name]=StringIO()
                else:
                    self.__Seq[name].write(line.strip())
        self.__Seq[name].seek(0)
        return toGetValue.keys()


    def readEMBL(self,emblStr):
        lines=emblStr.split("XX")

        try:
            accession=[x.strip()[2:].strip() for x in lines if x.strip()[:2]=='AC'][0]
            accession=accession.split(";")[0].split(",")[0]
            seq=[x for x in lines if x.strip()[:2]=='SQ']
            seq=seq[0]
            seqParts=seq.strip().split("\n")
            seq="".join(seqParts[1:])
            description="XX".join(lines[:-1]+["\n"+seqParts[0]])
            seqParts=None
        except IndexError:
            raise UnsupportedTypeException()

        seq=seq.translate(string.maketrans("",""),"/\n\r\t 0123456789")
        self.__Desc[accession]=description
        self.__Seq[accession]=StringIO(seq)

        return [accession]
        
        
    
    def addSequence(self, filename):
        "adds Sequences from file"
        try:
            # Slightly cludgy system to make it work with named pipes (FIFOs)
            # Can only read and open once and can't seek at all.
            File=open(filename,"r")
            fdata=File.read()
            File.close()
            try:
                GZFile=GzipFile(None,fileobj=StringIO(fdata))
                fileStr=GZFile.read()
            except (NameError,IOError):
                # NameError for missing GzipFile class
                # IOError for uncompressed data
                fileStr=fdata
            name=''
            if fileStr:
                if fileStr[0]=='>':
                    toGetValue=self.readFasta(fileStr)
                elif len(fileStr)>1 and fileStr[:2]=='ID':
                    toGetValue=self.readEMBL(fileStr)
                else:
                    raise UnsupportedTypeException()

            for name in toGetValue:
                self.__Seq[name].seek(0)
                self.__Seq[name]=(SingleSequence(self.__Seq[name]),0)


            print filename, "added" 
        except UnsupportedTypeException,e:
            e.args=(filename+"in not in supported format!",)
            raise e

    def __str__(self):
        """returns the names of the sequences"""
        return "\n".join(["%s : %s"%(name,self.describeShort(name)) for name in self.__Seq.keys()])
##        outs=""
##        for name in self.__Seq.keys():
##            outs+=name+"\n"
##        return outs
##        for name,seq in self.__Seq.items():
##            Ac,Cc,Gc,Tc=seq.count("A")+seq.count("a"),seq.count("C")+seq.count("c"),seq.count("G")+seq.count("g"),seq.count("T")+seq.count("t")
##            Nc=len(seq)-Ac-Tc-Cc-Gc
##            tot=float(len(seq))/100.0
##            outs+="%s\nLen A C G T Other= %d %4.2f %4.2f %4.2f %4.2f %4.2f\n\n"%(name,len(seq),Ac/tot,Cc/tot,Gc/tot,Tc/tot,Nc/tot)
        return outs


    def removeSequence(self, namePat):
        "Remove sequence by name"
        import re
        pat=re.compile(namePat)
        seqToRemove=[x for x in  self.__Seq.keys() if pat.match(x)]
        for name in seqToRemove:
            del self.__Seq[name]
        return seqToRemove
    
    def getNames(self):
        "Return a list of sequence names"
        return self.__Seq.keys()

    def keys(self):
        """Synonym for getNames()"""
        return self.getNames()
    
    def has_key(self,key):
        "True if have sequence called key"
        return self.__Seq.has_key(key)
    
    def __getitem__(self,key):
        """Returns the DNA sequence of sequence key.

        Returns actually a SigleSequence object"""
        return self.__Seq[key][0]


    def __len__(self):
        """Returns the number of loaded sequences"""
        return len(self.__Seq)

    def sequence(self,name):
        "Returns the DNA sequence of sequence name"
        return self[name]

    def describe(self,name):
        "Returns the description of the sequence"
        return self.__Desc[name].strip()

    def describeShort(self,name):
        "Return 1 line description of the sequence"
        return self.__Desc[name].split("\n")[0]
