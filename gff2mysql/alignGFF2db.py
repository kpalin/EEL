
### 
### $Log$
### Revision 1.1  2004/04/29 13:59:46  kpalin
### Initial revision
###
###


import MySQLdb
import sys,xreadlines
import re



class DBloader:
    def __init__(self):
        "Minimal constructor"
        self._modules=[]
        self._attrPat=re.compile(r"([A-Za-z][A-Za-z0-9_]+) (.*)")
        self._multiAlign=[]
        self.totLen=0
        self.coordReg=re.compile(r'(\w+)\|(\d+).(\d+)-(\d+)')
        self.db = MySQLdb.connect(host="localhost", user="bbu-palin", passwd="palindrome", db="bbu_palin")

    def getSeqId(self,species,seqName,dbName=None):
        cur=self.db.cursor()

        q="SELECT id FROM sequence WHERE species='%s' AND name='%s'"%(species,seqName)
        print q
        c=cur.execute(q)
        if c==0:
            cur.execute("INSERT INTO sequence (species,name) VALUES ('%s','%s')"%(species,seqName))
            c=cur.execute(q)
            print q
        seqs=cur.fetchall()
        print seqs
        assert(len(seqs)>=1)
        return int(seqs[0][0])



    def getTransFactorId(self,name,width):
        cur=self.db.cursor()

        q="SELECT id FROM tfactor WHERE name='%s' AND width=%d"%(name,width)
        print q
        c=cur.execute(q)
        if c==0:
            cur.execute("INSERT INTO tfactor (name,width) VALUES ('%s',%d)"%(name,width))
            c=cur.execute("SELECT LAST_INSERT_ID()")
        seqs=cur.fetchall()
        print seqs
        assert(len(seqs)>=1)
        return int(seqs[0][0])




    def _parseAttribs(self,attr):
        """Parses the attributes on GFF file line.

        Follows the attribute definition given in the standard
        from the Sanger center. Returns the mapping from attribute
        names to their values."""
        ret={}
        for a in attr.split(";"):
            try:
                attName,attVal=self._attrPat.match(a.strip()).groups()
                ret[attName]=attVal
            except AttributeError:
                pass
        return ret


    def parseSeqName(self,sname):
        m=self.coordReg.match(sname)
        if not m:
            raise "NO GOOD SEQUENCE NAME"
        vals={}
        vals["species"]=m.group(1)
        vals["name"]=m.group(2)
        vals["start"]=int(m.group(3))
        vals["end"]=int(m.group(4))
        return vals
    

    def addRegion(self,seqId,start,stop,cisId):
        cur=self.db.cursor()
        cur.execute("INSERT INTO region (seqID,beginPos,endPos,cisID) VALUES (%d,%d,%d,%d)"%(seqId,start,stop,cisId))
        cur.execute("SELECT LAST_INSERT_ID()")
        rid=cur.fetchall()
        assert(len(rid)==1)
        return int(rid[0][0])

    
    def addCisModule(self, score,s1,s1start,s1stop,s2,s2start,s2stop):

        seq1vals=self.parseSeqName(s1)
        seq1id=self.getSeqId(seq1vals["species"],seq1vals["name"])
        
        seq2vals=self.parseSeqName(s2)
        seq2id=self.getSeqId(seq2vals["species"],seq2vals["name"])

        cur=self.db.cursor()
        cur.execute("INSERT INTO cismodule (score) VALUES (%f)"%(score))
        cisId=int(cur.execute("SELECT LAST_INSERT_ID()").fetchall()[0][0])

        regionMap={}
        
        reg1=self.addRegion(seq1id,seq1vals["start"]+s1start,seq1vals["start"]+s1stop,cisId)
        regionMap[s1]=(reg1,seq1vals["start"])
        
        reg2=self.addRegion(seq2id,seq2vals["start"]+s2start,seq2vals["start"]+s2stop,cisId)
        regionMap[s2]=(reg2,seq2vals["start"])


        return regionMap
        
    def _parseGFFfile(self,fname):
        """Parses the pairwise alignment GFF file.

        This is used to clear up the addGFFfile() function."""
        fhandle=open(fname)

        currMod=None
        currModId=""
        cisModRows=0
        for line in xreadlines.xreadlines(fhandle):
            line=line.strip()
            if len(line)==0 or line[0]=='#':  ## Skip empty and comment lines.
                continue
            ## Parse a GFF line
            try:
                seq,src,feat,start,stop,score,strand,frame,attribs=line.split('#')[0].split("\t",8)
            except ValueError:
                print line
                raise
            attribs=self._parseAttribs(attribs)
            # Start of a new cis module
            if feat=='CisModule':
                try:
                    seq2=attribs["Target"].strip('"')
                except TypeError:
                    seq2=attribs["Target"].replace('"','')
                start2=int(attribs["Start"])
                stop2=int(attribs["End"])
                seqsHaveCol={}
                
                if cisModRows<1:
                    regionMap=self.addCisModule(float(score),seq,int(start),int(stop),\
                                                seq2,start2,stop2)
                    cisModRows+=1
            else:
                cisModRows=0

                if seqsHaveCol.has_key(seq):
                    seqsHaveCol[seq]=1
                    colId= UUSI
                cur=self.db.cursor()
                ins="INSERT INTO sites (pos,regID,colID) VALUES (%d,%d,%d)"%(regionMap[seq][1]+int(start),regionMap[seq][0],colId)
                cur.execute(ins)
                #assert(lineModId==currModId)
                #siteId=Site(feat,float(score),strand,int(stop)-int(start))
                #currMod.appendSite({seq:int(start)},siteId)
        if currMod:
            self._modules.append(currMod)

    def addGFFfile(self,fname):
        """Add a GFF file as input for the multiple alignment.

        The given GFF file should follow the format produced
        by mabs align command."""
        self._parseGFFfile(fname)
                


a=DBloader()
a.addGFFfile(sys.argv[1])
