
### 
### $Log$
### Revision 1.1  2004/04/29 14:01:14  kpalin
### Initial revisions
###
###


import MySQLdb
import sys,xreadlines
import re,gzip

HOSTNAME="localhost"
USERNAME="bbu-palin"
PASSWORD="palindrome"
DATABASE="bbu_palin"

SPECIES2DB={"HUMAN":"homo_sapiens_core_19_34a",
            "MOUSE":"mus_musculus_core_19_30"}



class DBloader:
    def __init__(self):
        "Minimal constructor"
        self._modules=[]
        self._attrPat=re.compile(r"([A-Za-z][A-Za-z0-9_]+) (.*)")
        self._multiAlign=[]
        self.totLen=0
        self.coordReg=re.compile(r'([^.]+).(\d+)-(\d+)')

        self.printIt=0
        self.resetConnection()

    def resetConnection(self):
        if hasattr(self,"db"):
            self.db.close()
        self.db = MySQLdb.connect(host=HOSTNAME, user=USERNAME, passwd=PASSWORD, db=DATABASE)
        self.seqCodeBuffer={}
        self.TFcodeBuffer={}
        self.addFileCount=0
        

    def getSeqId(self,species,seqName,dbName="default"):

        seqKey=(species,seqName,dbName)
        if not self.seqCodeBuffer.has_key(seqKey):
            cur=self.db.cursor()

            q="SELECT id FROM sequence WHERE species='%s' AND name='%s'"%(species,seqName)
            print q
            c=cur.execute(q)
            if c==0:
                cur.execute("INSERT INTO sequence (species,name,DBname) VALUES ('%s','%s','%s')"%(species,seqName,dbName))
                c=cur.execute("SELECT LAST_INSERT_ID()")
            seqs=cur.fetchall()
            print seqs
            self.printIt=1
            assert(len(seqs)>=1)

            self.seqCodeBuffer[seqKey]=int(seqs[0][0])
        return self.seqCodeBuffer[seqKey]


    def getTFid(self,tfName,width):
        tfKey=(tfName,width)
        if not self.TFcodeBuffer.has_key(tfKey):
            cur=self.db.cursor()

            q="SELECT id FROM tfactor WHERE name='%s' AND width=%d"%(tfName,width)

            c=cur.execute(q)
            if c==0:
                cur.execute("INSERT INTO tfactor (name,width) VALUES ('%s',%d)"%(tfName,width))
                c=cur.execute("SELECT LAST_INSERT_ID()")
            tfID=cur.fetchall()
            self.TFcodeBuffer[tfKey]=int(tfID[0][0])

        return self.TFcodeBuffer[tfKey]


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
        parts=sname.split("|")

        m=filter(lambda x:x,map(self.coordReg.match,parts))
        if len(m)==0:
            raise ("NO GOOD SEQUENCE NAME",sname)
        elif len(m)>1:
            print "Something fishy:",parts
            print "Assuming",m[0].groups()
        m=m[0]
            

        vals={}
        vals["species"]=parts[0]
        
        if "Crick" in parts:
            vals["strand"]="Crick"
        else:
            vals["strand"]="Watson"
            
        
        vals["name"]=m.group(1)
        vals["start"]=int(m.group(2))
        vals["end"]=int(m.group(3))
        #print vals["species"],vals["strand"],vals["name"],":",sname
        return vals
    

    def addRegion(self,seqId,start,stop,cisId):
        cur=self.db.cursor()
        cur.execute("INSERT INTO region (seqID,beginPos,endPos,cisID) VALUES (%d,%d,%d,%d)"%(seqId,start,stop,cisId))
        cur.execute("SELECT LAST_INSERT_ID()")
        rid=cur.fetchall()
        assert(len(rid)==1)
        return int(rid[0][0])

    
    def addCisModule(self, score,s1,s1start,s1stop,s2,s2start,s2stop):
        #print "New module",score

        seq1vals=self.parseSeqName(s1)
        seq1id=self.getSeqId(seq1vals["species"],seq1vals["name"],SPECIES2DB[seq1vals["species"]])
        if self.printIt:
            print s1
            self.printIt=0
        
        seq2vals=self.parseSeqName(s2)
        seq2id=self.getSeqId(seq2vals["species"],seq2vals["name"],SPECIES2DB[seq2vals["species"]])
        if self.printIt:
            print s2
            self.printIt=0

        cur=self.db.cursor()
        cur.execute("INSERT INTO cismodule (score) VALUES (%f)"%(score))
        cur.execute("SELECT LAST_INSERT_ID()")
        cisId=int(cur.fetchall()[0][0])

        regionMap={}
        
        reg1=self.addRegion(seq1id,seq1vals["start"]+s1start,seq1vals["start"]+s1stop,cisId)
        regionMap[s1]=(reg1,seq1vals)
        
        reg2=self.addRegion(seq2id,seq2vals["start"]+s2start,seq2vals["start"]+s2stop,cisId)
        regionMap[s2]=(reg2,seq2vals)


        return regionMap


    def makeNewColumn(self,tfName,width,scoreDelta):
        cur=self.db.cursor()

        tfID=self.getTFid(tfName,width)
        ins="INSERT INTO alnCols (tfID,scoredelta) VALUES (%d,%g)"%(tfID,scoreDelta)
        cur.execute(ins)
        cur.execute("SELECT LAST_INSERT_ID()")
        rid=cur.fetchall()
        colId=int(rid[0][0])
        return colId
        
    
    def _parseGFFfile(self,fhandle):
        """Parses the pairwise alignment GFF file.

        This is used to clear up the addGFFfile() function."""

        currMod=None
        currModId=""
        cisModRows=0

        # Book-keeping for columns
        seqsHaveCol={}

        
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
                
                if cisModRows<1:
                    regionMap=self.addCisModule(float(score),seq,int(start),int(stop),\
                                                seq2,start2,stop2)
                    cisModRows+=1

                # Must create the first column.
                seqsHaveCol[seq]=1
            else:
                cisModRows=0

                if seqsHaveCol.has_key(seq):
                    # New column
                    seqsHaveCol={}
                    colId=self.makeNewColumn(feat,int(stop)-int(start),0.0)
                seqsHaveCol[seq]=1

                
                cur=self.db.cursor()

                # To get the strands correctly:
                if regionMap[seq][1]["strand"]=="Crick":
                    if strand=='-':
                        strand='+'
                    elif strand=='+':
                        strand='-'
                    pos=regionMap[seq][1]["end"]-int(stop)
                else:
                    pos=regionMap[seq][1]["start"]+int(start)

                    
                ins="INSERT INTO sites (pos,regID,colID,strand) VALUES (%d,%d,%d,'%s')"%(pos,regionMap[seq][0],colId,strand)
                cur.execute(ins)
                #assert(lineModId==currModId)
                #siteId=Site(feat,float(score),strand,int(stop)-int(start))
                #currMod.appendSite({seq:int(start)},siteId)

        print "Done"

            
    def addGFFfile(self,fname):
        """Add a GFF file as input for the multiple alignment.

        The given GFF file should follow the format produced
        by mabs align command."""
        self.addFileCount+=1
        if self.addFileCount>10:
            self.resetConnection()
            self.addFileCount=0
            
        try:
            fhandle=gzip.GzipFile(fname)
            self._parseGFFfile(fhandle)
        except IOError:
            fhandle=open(fname)
            self._parseGFFfile(fhandle)
        print "Added",fname
                


a=DBloader()
import glob

for fname in glob.glob(sys.argv[1]):
    a.addGFFfile(fname)
