
### 
### $Log$
### Revision 1.2  2004/05/04 10:26:43  kpalin
### Toimii jotenkin.
###
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

resetConnectionAfter=100

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
        self.sliceCodes={}
        

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
            cur.close()
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
            cur.close()
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
        cur.close()
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
        vals["name"]=m.group(1)
        vals["start"]=int(m.group(2))
        vals["end"]=int(m.group(3))


        vals["species"]=parts[0]
        
        if "Crick" in parts:
            vals["strand"]="Crick"
        else:
            vals["strand"]="Watson"
            
        geneid=[x for x in parts if x[:3]=='ENS']
        if len(geneid)==1:
            vals["gene"]=geneid[0]
        elif len(geneid)>1:
            print "Funny",geneid
        
        #print vals["species"],vals["strand"],vals["name"],":",sname
        return vals
    

    def addSlice(self,cisID,regionMap,fname):
        cur=self.db.cursor()
        for seq in regionMap.keys():

            if self.sliceCodes.has_key((fname,seq)):
                sliceID=self.sliceCodes[(fname,seq)]
            else:
                c=cur.execute("SELECT id FROM slice WHERE fileName='%s' and seqName='%s'"%(fname,seq))
                if c==1:
                    sliceID=int(cur.fetchall()[0][0])
                elif c==0:
                    if regionMap[seq][1]["strand"]=="Crick":
                        revComp=1
                    else:
                        revComp=0
                    start=regionMap[seq][1]["start"]
                    stop=regionMap[seq][1]["end"]

                    geneName=""
                    if regionMap[seq][1].has_key("gene"):
                        geneName=regionMap[seq][1]["gene"]
                    cur.execute("INSERT INTO slice (beginPos,endPos,fileName,seqName,geneName,revComplement) VALUES (%d,%d,'%s','%s','%s',%d)"%(start,stop,fname,seq,geneName,revComp))
                    cur.execute("SELECT LAST_INSERT_ID()")
                    sliceID=int(cur.fetchall()[0][0])
                else:
                    print "VERY WEIRD:"
                    print cur.fetchall()
                    
                self.sliceCodes[(fname,seq)]=sliceID

            cur.execute("INSERT INTO sliceCisLink (sliceID,cisID) VALUES (%d,%d)"%(sliceID,cisID))
        cur.close()

    def addRegion(self,seqId,start,stop,cisId):
        cur=self.db.cursor()
        cur.execute("INSERT INTO region (seqID,beginPos,endPos,cisID) VALUES (%d,%d,%d,%d)"%(seqId,start,stop,cisId))
        cur.execute("SELECT LAST_INSERT_ID()")
        rid=cur.fetchall()
        assert(len(rid)==1)
        cur.close()
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


        self.currentCisID=cisId

        regionMap={}
        s1start,s1stop,strand=self.mapWithinRegion(seq1vals,s1start,s1stop)
        assert(s1start<s1stop)

        reg1=self.addRegion(seq1id,seq1vals["start"]+s1start,seq1vals["start"]+s1stop,cisId)
        regionMap[s1]=(reg1,seq1vals)
        
        s2start,s2stop,strand=self.mapWithinRegion(seq1vals,s2start,s2stop)
        assert(s2start<s2stop)

        reg2=self.addRegion(seq2id,seq2vals["start"]+s2start,seq2vals["start"]+s2stop,cisId)
        regionMap[s2]=(reg2,seq2vals)

        cur.close()
        return regionMap


    def makeNewColumn(self,tfName,width,scoreDelta):
        cur=self.db.cursor()

        tfID=self.getTFid(tfName,width)
        ins="INSERT INTO alnCols (tfID,scoredelta) VALUES (%d,%g)"%(tfID,scoreDelta)
        cur.execute(ins)
        cur.execute("SELECT LAST_INSERT_ID()")
        rid=cur.fetchall()
        colId=int(rid[0][0])
        cur.close()
        return colId
        

    def mapWithinRegion(self,regionMap,start,stop,strand='.'):
        # To get the strands correctly:
        if regionMap["strand"]=="Crick":
            if strand=='-':
                strand='+'
            elif strand=='+':
                strand='-'

            stop=regionMap["end"]-int(start)
            if int(stop)>=0:
                start=regionMap["end"]-int(stop)
            else:
                start=regionMap["start"]-int(stop)-1
        else:
            start=regionMap["start"]+int(start)
            if int(stop)>=0:
                stop=regionMap["start"]+int(stop)
            else:
                stop=regionMap["end"]+int(stop)+1  # stop==-1 <=> last

        return (start,stop,strand)
    
    def _parseGFFfile(self,fhandle,fileName):
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

                    self.addSlice(self.currentCisID,regionMap,fileName)
                
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

                start,stop,strand=self.mapWithinRegion(regionMap[seq][1],start,stop,strand)
                assert(start<stop)

                ins="INSERT INTO sites (pos,regID,colID,strand) VALUES (%d,%d,%d,'%s')"%(start,regionMap[seq][0],colId,strand)
                cur.execute(ins)
                #assert(lineModId==currModId)
                #siteId=Site(feat,float(score),strand,int(stop)-int(start))
                #currMod.appendSite({seq:int(start)},siteId)
        cur.close()
        print "Done"

            
    def addGFFfile(self,fname):
        """Add a GFF file as input for the multiple alignment.

        The given GFF file should follow the format produced
        by mabs align command."""
        self.addFileCount+=1
        if self.addFileCount>resetConnectionAfter:
            self.resetConnection()
            self.addFileCount=0
            
        try:
            fhandle=gzip.GzipFile(fname)
            self._parseGFFfile(fhandle,fname)
        except IOError:
            fhandle=open(fname)
            self._parseGFFfile(fhandle,fname)
        print "Added",fname
                



def main():
    a=DBloader()
    import glob
    flist=glob.glob(sys.argv[1])
    for find in range(len(flist)):
        print "Adding %d/%d"%(find+1,len(flist))
        fname=flist[find]
        a.addGFFfile(fname)


if 0:
    import hotshot
    import hotshot.stats

    profiler= hotshot.Profile("logfile.dat")
    profiler.run("main()")
    profiler.close()
    stats=hotshot.stats.load("logfile.dat")

    stats.strip_dirs().sort_stats("cumulative").print_stats("alignGFF")
    
else:
    main()
