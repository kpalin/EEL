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
    pass

class Sequences:
    "represents DNA-Sequences"
    def __init__(self,filename=0):
        "reads seqences from file"
        self.__Seq={}
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
                    name=string.split(line[1:].strip())[0]
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
            seq="".join(seq.strip().split("\n")[1:])
        except IndexError:
            raise UnsupportedTypeException()

        seq=seq.translate(string.maketrans("",""),"/\n\r\t 0123456789")
        self.__Seq[accession]=StringIO(seq)

        return [accession]
        
        
    
    def addSequence(self, filename):
        "adds Sequences from file"
        try:
            try:
                File=GzipFile(filename,'r')
                File.read(1)
                File.seek(0)
            except (NameError,IOError):
                File=open(filename,"r")
            name=''
            fileStr=File.read()
            File.close()
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
        except IOError, (errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        except UnsupportedTypeException:
            print filename,"in not in supported format!"

    def __str__(self):
        """returns the names of the sequences"""
        outs=""
        for name in self.__Seq.keys():
            outs+=name+"\n"
        return outs
        for name,seq in self.__Seq.items():
            Ac,Cc,Gc,Tc=seq.count("A")+seq.count("a"),seq.count("C")+seq.count("c"),seq.count("G")+seq.count("g"),seq.count("T")+seq.count("t")
            Nc=len(seq)-Ac-Tc-Cc-Gc
            tot=float(len(seq))/100.0
            outs+="%s\nLen A C G T Other= %d %4.2f %4.2f %4.2f %4.2f %4.2f\n\n"%(name,len(seq),Ac/tot,Cc/tot,Gc/tot,Tc/tot,Nc/tot)
        return outs


    def removeSequence(self, name):
        "Remove sequence by name"
        if self.__Seq.has_key(name):
            del self.__Seq[name]
        
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
