import string
from cStringIO import StringIO
from xreadlines import xreadlines

from gzip import GzipFile


class SingleSequence:
    def __init__(self,file):
        self.file=file

        self.start=self.file.tell()
        self.file.seek(0,2)

        self.end=self.file.tell()
        self.file.seek(self.start,0)

        self.len=self.end-self.start
        
    def __getattr__(self,attr):
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
        c=self.file.read(i-j)
        self.file.seek(self.start)

        return c
    
    def __len__(self):
        return self.len


class Sequences:
    "represents DNA-Sequences"
    def __init__(self,filename=0):
        "reads seqences from file"
        self.__Seq={}
        if filename:
            self.addSequence(filename)



    def resetSeqPositions(self):
        for seq,pos in self.__Seq.values():
            seq.seek(pos)

    def addSingleSequence(self,filename):
        """Adds a single sequence from fasta formated and gziped file.
        Useful for huge sequences, e.g. whole chromosomes"""
        try:
            File=GzipFile(filename,'r')
            line="ALUSSA AINA HANKALAA"
            try:
                line=File.readline()
            except IOError,e:
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


        
    def addSequence(self, filename):
        "adds Sequences from file"
        try:
            File=open(filename,'r')
            name=''
            #for line in string.split(File.read(),"\n"):
            toGetValue={}
            for line in xreadlines(File):
                if line and line[0]=='>':
                    name=string.split(line[1:].strip())[0]
                    toGetValue[name]=1
                    #name=line[1:].strip()
                    self.__Seq[name]=StringIO()
                else:
                    self.__Seq[name].write(line.strip())
            File.close()
            for name in toGetValue.keys():
                self.__Seq[name].seek(0)
                self.__Seq[name]=(SingleSequence(self.__Seq[name]),0)


            print filename, "added" 
        except IOError, (errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        except KeyError:
            print filename,"in not in FASTA format!"

    def __str__(self):
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
        if self.__Seq.has_key(name):
            del self.__Seq[name]
        
    def getNames(self):
        return self.__Seq.keys()

    def has_key(self,key):
        return self.__Seq.has_key(key)
    
    def __getitem__(self,key):
        return self.__Seq[key][0]

    def sequence(self,name):
        return self.__Seq[name][0]
