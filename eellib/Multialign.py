"""Python interface and logic for multiple alignment."""

#
# $Log$
# Revision 1.8  2004/04/08 13:03:12  kpalin
# Output and stupid Greedy multiple alignment works.
#
# Revision 1.7  2004/04/06 12:18:34  kpalin
# Looks a lot like the Prohaska way doesn't work. Trying something
# nasty.
#
# Revision 1.6  2004/03/29 12:28:41  kpalin
# Should work for three species but not more. Trying to fix
# the more multiple stuff.
#
# Revision 1.5  2004/03/03 09:16:21  kpalin
# Functional version with reasonably good integration to mabs interface
#
# Revision 1.4  2004/03/02 10:45:56  kpalin
# Working version. It doesn't give proper output yet and it's not
# integrated to the mabs interface but it mostly works.
#
# Revision 1.3  2004/03/01 10:19:18  kpalin
# Starting to implement an other aproach by Prohaska et.al.
#
# Revision 1.2  2004/02/23 12:23:52  kpalin
# Updates for per gene orthologous runs. Maybe litle multiple alignment.
#
# Revision 1.1  2004/02/19 13:17:51  kpalin
# Documentation for the main parts of multiple alignment
# done. Now we only need the code.
#
#


import xreadlines,re,traceback,pdb,math
from StringIO import StringIO

class NotComparable(Exception):
    def __init__(self,f,s):
        self.first=f
        self.second=s

    def __str__(self):
        return "Can't compare %s <=> %s"%(str(self.first),str(self.second))

class PairModule:
    """A pairwise cis module from a GFF file.

    Use as a representation of the pairwise alignment cis
    module from the mabs alignment GFF files. This class
    has several important functions for adding pairs to form
    multiple alignments."""

    def __init__(self,cisCode,score,sname1,spos1,epos1,sname2,spos2,epos2):
        """Constructor defining the sequence name and position
        of the cis module.

        cisCode is the unique code for this cis module. score
        is the total pairwise alignment score.
        The parameter sname1 is the name for the first sequence,
        spos1 is the start and epos is the end position on the
        first sequence DNA. The parameters sname2, spos2 and epos2
        are the same for the second sequence."""
        self.id=cisCode
        self.score=score

        self.sNames=[sname1,sname2]
        self.sNames.sort()
        # For greedy alignment.
        self.used=-1
        self.included=0
        
        self.DNApos={}
        self.DNApos[sname1]=(spos1,epos1)
        self.DNApos[sname2]=(spos2,epos2)
        self.sitePos={sname1:[],sname2:[]}
        self.sites=[]

    def __cmp__(self,other):
        """Comparison between two PairModule:s

        The PairModules are partially ordered along their alphabetically
        first common sequence."""
        haveCommon=0
        for common in self.sNames:
            if common in other.sNames:
                haveCommon=1
                break
        #try:
        #    assert(haveCommon==1)
        #except AssertionError:
        #    raise NotComparable(self,other)

        #try:
        r=0
        if haveCommon:
            r=cmp(self.DNApos[common],other.DNApos[common])
        #except KeyError:
        #    print common
        #    print self.DNApos
        #    print other.DNApos
        #    raise
        if r==0 or not haveCommon:
            r=cmp(self.sNames,other.sNames)

        return r

    def intersects(self,other):
        """Checks if there is intersection between the two PairModules.

        Return true if the modules share a common TF binding site,
        false otherwise."""
        if not intersection(self,other):
            return 0
        return 1

    def intersection(self,other,common=None):
        """Returns the intersection of two PairModules.

        Return value is the list of pairs of overlaping sites."""
        #print "Intersecting",self.sNames,other.sNames
        if not common:
            for common in self.sNames:
                if other.DNApos.has_key(common):
                    break
        try:
            # No overlap if self ends before beginning of other
            # or starts after the end of the other
            if self.DNApos[common][1]<other.DNApos[common][0] or \
               self.DNApos[common][0]>other.DNApos[common][1]:
                return 0
        except KeyError:
            return 0

        i,j=0,0
        inters=[]
        while i<len(self.sites) and j<len(other.sites):
            if self.sitePos[common][i]<other.sitePos[common][j]:
                i+=1
            elif self.sitePos[common][i]>other.sitePos[common][j]:
                j+=1
            elif self.sites[i].match(other.sites[j]):
                inters.append((i,j))
                i+=1
                j+=1
            else:
                i+=1
                j+=1
        if len(inters)==0:
            return 0
        return inters
                
    def makeTrace(self,other):
        """Make a trace. That is a new pairModule subset of self where other
        is projected."""


        sname1,sname2=self.sNames
        trace=None
        for common in self.sNames:
            try:
                # No overlap if self ends before beginning of other
                # or starts after the end of the other
                if self.DNApos[common][1]<other.DNApos[common][0] or \
                   self.DNApos[common][0]>other.DNApos[common][1]:
                    continue
            except KeyError:
                continue
            if not trace:
                trace=PairModule(None,None,sname1,0,1,sname2,0,1)
            i,j=0,0
            while i<len(self.sites) and j<len(other.sites):
                if self.sitePos[common][i]<other.sitePos[common][j]:
                    i+=1
                elif self.sitePos[common][i]>other.sitePos[common][j]:
                    j+=1
                elif self.sites[i].match(other.sites[j]):
                    trace.sites.append((self.sitePos[sname1][i],self.sites[i]))
                    trace.sitePos[sname1].append(self.sitePos[sname1][i])
                    trace.sitePos[sname2].append(self.sitePos[sname2][i])
                    i+=1
                    j+=1
                else:
                    i+=1
                    j+=1

        if trace and len(trace.sites)>0:
            trace.sites.sort()
            trace.sitePos[sname1].sort()
            trace.sitePos[sname2].sort()
            trace.sites=[x[1] for x in trace.sites]
            trace.DNApos[sname1]=(trace.sitePos[sname1][0],trace.sitePos[sname1][-1]+trace.sites[-1].width)
            trace.DNApos[sname2]=(trace.sitePos[sname2][0],trace.sitePos[sname2][-1]+trace.sites[-1].width)

            trace.oldPos=self.DNApos
            trace.score=0.0
            trace.id="trace"
        else:
            trace=None
            
        return trace
            
            
    def isSubOrNoCommon(self,other):
        """Return true if other is subalignment of self or they don't share
        common sequence"""
        cc=0
        fullMatch=0
        ennen=-1
        for common in self.sNames:
            try:
                # No overlap if self ends before beginning of other
                # or starts after the end of the other
                if self.DNApos[common][1]<other.DNApos[common][0] or \
                   self.DNApos[common][0]>other.DNApos[common][1]:
                    return 0
            except KeyError:
                continue
            cc+=1
            i,j=0,0
            while i<len(self.sites) and j<len(other.sites):
                if self.sitePos[common][i]<other.sitePos[common][j]:
                    return 0
                    i+=1
                elif self.sitePos[common][i]>other.sitePos[common][j]:
                    j+=1
                elif self.sites[i].match(other.sites[j]):
                    i+=1
                    j+=1
                else:
                    return 0
            if i==len(self.sites):
                fullMatch+=1

        if fullMatch>0:
            return 1
        else:
            print self.DNApos.keys(),other.DNApos.keys(),common
            print "Not full"
            return 0

    def isSubOrNoCommon2(self,other):
        """Return true if other is subalignment of self or they don't share
        common sequence"""

        for common in [x for x in self.sNames if x in other.sNames]:
            i,j=0,0
            while i<len(self.sites) and j<len(other.sites):
                if self.sitePos[common][i]<other.sitePos[common][j]:
                    return 0
                elif self.sitePos[common][i]>other.sitePos[common][j]:
                    j+=1
                elif self.sites[i].match(other.sites[j]):
                    i+=1
                    j+=1
                else:
                    return 0
            if i<len(self.sites):
                return 0
        return 1

    
    def appendSite(self,sPos,site):
        """Append a new pair of sites belonging to this cis module.

        The parameter sPos is a map from sequence name to tuple of sites
        DNA coordinate and match score. site is the TF binding site
        identification."""
        for seq,pos in sPos.items():
            if seq==self.sNames[0]:
                self.sitePos[seq].append(pos)
                self.sites.append(site)
            else:
                self.sitePos[seq].append(pos)                
        pass

    def _finalize(self):
        assert(None)
        for seq in self.sNames:
            self.sitePos[seq].sort()


    def __str__(self):
        #self._finalize()
        
        strbuf=StringIO()
        
        cmFormat='%s\tmalign\tCisModule\t%d\t%d\t%g\t.\t.\tTarget "%s";\tStart %d;\tEnd %d;\tStrand .;\tFrame .;\tCM %s;\n'

        seq1,seq2=self.sNames
        strbuf.write(cmFormat%(seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.score,seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.id))
        strbuf.write(cmFormat%(seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.score,seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.id))

        return strbuf.getvalue()

    def toGFF(self):
        """Returns GFF representation of the pairModule."""
        #self._finalize()
        
        strbuf=StringIO()
        cmFormat='%s\tmalign\tCisModule\t%d\t%d\t%g\t.\t.\tTarget "%s";\tStart %d;\tEnd %d;\tStrand .;\tFrame .;\tCM %s;\n'
        seq1,seq2=self.sNames
        strbuf.write(cmFormat%(seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.score,seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.id))
        strbuf.write(cmFormat%(seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.score,seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.id))

        siteFormat='%s\tmalign\t%s\t%d\t%d\t%g\t.\t.\tCM %s;\n'
        for i in range(len(self.sitePos[seq1])):
            pos1=self.sitePos[seq1][i]
            site=self.sites[i]
            strbuf.write(siteFormat%(seq1,site.name,pos1,pos1+site.width,site.score,self.id))
            pos2=self.sitePos[seq2][i]
            strbuf.write(siteFormat%(seq2,site.name,pos2,pos2+site.width,site.score,self.id))
        return strbuf.getvalue()
                         


class MultiOutModule:
    """Class to assist the output of multiple alignments"""
    def __init__(self,aligns,sites,colCounts,id,alnLimit=2,score=0.0,pairCount=0):
        self.score=score
        self.align=aligns
        self.sites=sites
        self.colCount=colCounts
        self.id=id
        self.n=len(self.sites)
        self.alnLimit=alnLimit
        self.pairCount=pairCount

        self.sb=None
        self.seqs=self.align.keys()
        self.seqContext=10
        self.seqCoords={}
        
        self.findSeqCoords()

    def setAlnLimit(self,alnLimit=0):
        self.sb=None
        self.alnLimit=alnLimit
        
    def findSeqCoords(self):
        self.seqCoords={}
        for seq,posses in self.align.items():
            for i in posses:
                if i>=0:
                    start=i
                    break
            for j in range(len(posses)-1,0,-1):
                if posses[j]>=0:
                    stop=posses[j]+self.sites[j].width
                    break
            self.seqCoords[seq]=(start,stop)
       
    def __len__(self):
        return len(self.sites)

    def writeSeqNames(self):
        self.sb.write("Sequences:\n")
        for seq in self.seqs:
            self.sb.write("%s\n"%(seq))
        self.sb.write("\n")

    def writePosAln(self):
        self.sb.write("Alignment Score: %g\n"%(self.score))
        if self.pairCount>0:
            self.sb.write("Used pairwise alignments: %d\n"%(self.pairCount))
        self.sb.write(" "*26+"  ".join(["%-10s"%(x[:10]) for x in self.seqs]))
        for i in range(self.n):
            #alnCount=(1+math.sqrt(1+8*self.colCount[i]))/2
            alnCount=self.colCount[i]

            if alnCount<self.alnLimit:
                continue
            self.sb.write("\n%-20s %s %2d "%(self.sites[i].name[-20:],self.sites[i].strand,int(alnCount)))
            self.sb.write("  ".join(["%10d"%(self.align[x][i]) for x in self.seqs]))
        self.sb.write("\n")

    def _readyToOutput(self):
        if not self.sb:
            self.sb=StringIO()
            self.writeSeqNames()
            self.writePosAln()
            
    def __str__(self):
        self._readyToOutput()
        return self.sb.getvalue()


    def writeAlnSeq(self,linelen=60):
        pos={}
        aln={}
        self.sb.write("\n")
        for (seq,i) in zip(self.seqs,range(len(self.seqs))):
            pos[seq]=self.seqCoords[seq][0]-self.seqContext
            aln[seq]=self.seqAln[seq].getvalue()
            self.sb.write("Sequence %d : %s\n"%(i,seq))

        alnLen=len(aln.values()[0])
        self.sb.write("\n")
        for p in range(0,alnLen,linelen):
            for seq in aln.keys():
                try:
                    alnLine=aln[seq][p:p+linelen]
                except KeyError:
                    print seq,self.seqs,aln.keys()
                    print self.sb.getvalue()
                    raise
                self.sb.write("%7d : %s\n"%(pos[seq],alnLine))
                pos[seq]=pos[seq]+len(alnLine)-alnLine.count("-")-alnLine.count(" ")
            self.sb.write("\n")
            
                
            
    def strAln(self,fullSeqs):
        """Format sequences with respect to the multiple alignment"""
        self._readyToOutput()
        prevPos={}
        myAlnLen={}
        self.seqAln={}
        times={}

        alnList=["."]
        
        for seq in self.seqs:
            prevPos[seq]=self.seqCoords[seq][0]-self.seqContext
            self.seqAln[seq]=StringIO()
            myAlnLen[seq]=0
            times[seq]=0

        alnLen=0
        for p in range(self.n):
            newAlnLen=alnLen
            #print p
            for seq in self.seqs:
                nextPos=self.align[seq][p]
                
                # Do nothing if this site is missing from this sequence
                if nextPos<0:
                    continue 

                # Otherwise,
                newAlnLen=max(newAlnLen,myAlnLen[seq]+nextPos-prevPos[seq])
            for seq in self.seqs:
                nextPos=self.align[seq][p]

                # Do nothing if this site is missing from this sequence
                if nextPos<0:
                    continue 

                # Otherwise,
                myInterval=nextPos-prevPos[seq]
                alnInterval=newAlnLen-myAlnLen[seq]

                if myInterval>0:
                    intSeq=fullSeqs[seq][prevPos[seq]:nextPos].lower()
                else:
                    intSeq=""

                if alnInterval>0:
                    if times[seq]>0:
                        # Justify the first interval (context) to the rigth
                        # others to the left.
                        apu=intSeq[:myInterval/2].ljust(alnInterval/2)+intSeq[myInterval/2:].rjust(int(math.ceil(alnInterval/2.0)))
                        intSeq=intSeq.ljust(alnInterval).replace(' ','-')
                        assert(len(apu)==len(intSeq)==alnInterval)
                        
                        intSeq=apu.replace(' ','-')
                    else:
                        intSeq=intSeq.rjust(alnInterval)
                    times[seq]+=1
                    self.seqAln[seq].write(intSeq)
                    myAlnLen[seq]+=len(intSeq)

                # Condition for site p not being contained in larger site starting before it
                if myInterval > -self.sites[p].width:
                    siteSeq=fullSeqs[seq][nextPos-min(0,myInterval):nextPos+self.sites[p].width].upper()
                    self.seqAln[seq].write(siteSeq)
                    myAlnLen[seq]+=len(siteSeq)

                    prevPos[seq]=nextPos+len(siteSeq)
                    #print seq[:10],self.seqAln[seq].getvalue()

            # Building alignment guide
            alnList.extend(["."]*(newAlnLen-len(alnList)+1))
            try:
                alnList[newAlnLen]=str(self.colCount[p])
            except Exception:
                pdb.set_trace()


        for seq in self.seqs:
            self.seqAln[seq].write(fullSeqs[seq][prevPos[seq]:prevPos[seq]+self.seqContext].lower())
        #raise "EHEI"

        self.seqAln["ALIGN"]=StringIO("".join(alnList)+"."*(max(myAlnLen.values())-len(alnList)+self.seqContext))
        self.seqCoords["ALIGN"]=(self.seqContext,max(myAlnLen.values()))

        self.writeAlnSeq()
                     
    def strAlnOld(self,fullSeqs):
        prevPos={}
        self.seqAln={}
        for seq in self.seqs:
            prevPos[seq]=self.seqCoords[seq][0]-self.seqContext
            self.seqAln[seq]=StringIO()

        #Give 10 bp of context before alignments.
        intervals=[self.seqContext]
        for p in range(1,self.n):
            intervals.append(-100)
            for seq in self.seqs:
                try:
                    if self.align[seq][p]>=0 and self.align[seq][p-1]>=0:
                        #                        if self.align[seq][p-1]>=0:
                        #                            thisInt=self.align[seq][p]-(self.sites[p-1].width+self.align[seq][p-1])
                        #                        else:
                        thisInt=-1e99
                        for j in range(p-1,-1,-1):
                            if self.align[seq][j]>=0:
                                thisInt=self.align[seq][p]-(self.sites[j].width+self.align[seq][j])
                                break
                        if thisInt<=0 or thisInt>1000:
                            pdb.set_trace()
                        intervals[-1]=max(intervals[-1],thisInt)
                except TypeError:
                    raise
                    pass

        print intervals
        # Align the sites. The non-site positions are left justified.
        for p in range(self.n):
            for seq in self.seqs:
                nextPos=self.align[seq][p]
                if nextPos<0:
                    continue
                #print "intervals[%d]=%d, nextPos-prevPos[seq]=%d"%(p,intervals[p],nextPos-prevPos[seq])
                thisInt=max(intervals[p],nextPos-prevPos[seq])
                intSeq=fullSeqs[seq][prevPos[seq]:nextPos].lower().ljust(thisInt).replace(" ","-")
                #print "%d. interval for %s: %d-%d : %d : %s"%(p,seq,prevPos[seq],nextPos,thisInt,intSeq)
                self.seqAln[seq].write(intSeq)
                self.seqAln[seq].write(fullSeqs[seq][nextPos:(nextPos+self.sites[p].width)].upper())
                prevPos[seq]=nextPos+self.sites[p].width

        # Last 10bp of context
        for seq in self.seqs:
            self.seqAln[seq].write(fullSeqs[seq][prevPos[seq]:prevPos[seq]+self.seqContext].lower())
            #print self.seqAln[seq].getvalue()
            #self.sb.write("\n")
        self.writeAlnSeq()

        


class MultiModule:
    """Class that forms to be the cis module in multiple sequences.

    Input is a long list of pair modules."""

    def __init__(self):
        """Minimal constructor.

        Most work is done while adding the pairwise alignments."""
        self._TFseqs={}
        self.sites=[]
        self._pairsBySeq={}
        self._pairs=[]

        # Adjacency matrix for overlap graph
        self._overlaps={}
        self._overlapsM={}
        pass


    def uniteButLast(self,other):
        """Unite self with other, excluding the last pair in other."""

        offset=len(self._pairs)
        last=len(other._pairs)-1

        self._pairs[-1].code=offset-1
        #print "I  ",[x.code for x in self._pairs],[x.code for x in other._pairs]
        #other._pairs[last].code=-1
        for op in other._pairs[:-1]:
            op.code=op.code+offset
            self._pairs.append(op)
        #self._pairs.extend(other._pairs[:-1])
        #print "II ",[x.code for x in self._pairs],[x.code for x in other._pairs[:-1]]
        newIndForLast=self._pairs.index(other._pairs[last])
        for (s,t),v in other._overlaps.items():
            
            if s==last:
                ns=newIndForLast
            else:
                ns=s+offset
            if t==last:
                nt=newIndForLast
            else:
                nt=t+offset
            if (ns<nt and s>t ) or (ns>nt and s<t):
                v=other._overlaps[(t,s)]
                if type(v)==type([]):
                    v=[(y,x) for (x,y) in v]
                
            self._overlaps[(ns,nt)]=v
            try:
                if not self._overlapsM.has_key(ns):
                    self._overlapsM[ns]={}
                self._overlapsM[ns][nt]=other._overlapsM[s][t]
            except (KeyError,TypeError):
                pass

        for name in other._pairsBySeq.keys():
            if not self._pairsBySeq.has_key(name):
                self._pairsBySeq[name]=[]
            self._pairsBySeq[name].extend([(pair,n+offset) for (pair,n) in filter(lambda x:x[1]!=last,other._pairsBySeq[name])])
        
        #print "unite",[x.code for x in self._pairs]
 
        
    def __str__(self):
        return str(self._pairsBySeq.items())


        

        
    def addPair(self,pair):
        """Adds a PairModule to this multiple alignment.

        The multiple alignment is formed by aligning all pairs
        of the sequences in the alignment with each other. These
        pairwise alignments are done before the multiple alignment."""


        n=len(self._pairs)
        # The first pairwise alignment.
        if n>0:

            self._overlapsM[n]={}
            
            # Do we have a common sequence?
            haveKey=0
            intersects=0
            for common in pair.sNames:
                if self._pairsBySeq.has_key(common):
                    # Do we have an overlapping pairwise alignment?
                    ## Bin search??
                    for opair,code in self._pairsBySeq[common]:
                        isects=opair.intersection(pair)
                        if isects:
                            intersects=1
                            self._overlaps[(code,n)]=isects
                            self._overlaps[(n,code)]=common
                            if not self._overlapsM[n].has_key(code):
                                self._overlapsM[n][code]=[common]
                            else:
                                self._overlapsM[n][code].append(common)
                            if not hasattr(pair,"overlaps"):
                                pair.overlaps=[code]
                            else:
                                pair.overlaps.append(code)
                            if not hasattr(opair,"overlaps"):
                                opair.overlaps=[n]
                            else:
                                opair.overlaps.append(n)
                                
            if not intersects:
                return 0

        for name in pair.sNames:
            if not self._pairsBySeq.has_key(name):
                self._pairsBySeq[name]=[]
            self._pairsBySeq[name].append((pair,n))

        pair.code=n
        assert(len(self._pairs)==n)
        self._pairs.append(pair)
        #print "add",[x.code for x in self._pairs],n
        assert(len(self._pairs)==(n+1))
        
        return 1                


    def _intersectingSitePairs(self,isect1,isect2):
        """Return the site indexes that can be aligned.

        Input is list of (a,b) and (b,c) where a,b,c are integers.
        The output is (a,c) such that b:s match."""
        # Find the intersecting sites.
        n,m=0,0
        isect=[]
        while n<len(isect1) and m<len(isect2):
            if isect1[n][1]<isect2[m][0]:
                n+=1
            elif isect1[n][1]>isect2[m][0]:
                m+=1
            else:
                assert(isect1[n][1]==isect2[m][0])
                npair=(isect1[n][0],isect2[m][1])
                isect.append(npair)
                n,m=n+1,m+1

        #print "i1",isect1
        #print "i2",isect2
        #print "i ",isect
        return isect

    def extendTransitivePath(self,i,j,k):
        """Extends the transitive path form i to j by k.

        This method is used in the inner most loop of the Floyd--Warshall
        algorithm in method TransitiveClosure(). The trick is to keep track
        of the exact mappings between aligned sites. Just tracking
        overlaps/doesn't overlap is not enough"""

        if i==j:
            if not self._TC.has_key((i,j)):
                l=len(self._pairs[i].sites)
                self._TC[(i,j)]=zip(range(l),range(l))
            return
        
        try:
            isectIJ=self._TC[(min(i,j),max(i,j))]
        except KeyError:
            isectIJ=[]
        isectIK=self._TC[(min(i,k),max(i,k))]
        isectKJ=self._TC[(min(j,k),max(j,k))]
        isect=self._intersectingSitePairs(isectIK,isectKJ)
        #print isect
        def _pairsAreSorted(l):
            r=1
            for i in range(1,len(l)):
                r=r and l[i-1][0]<l[i][0] and l[i-1][1]<l[i][1]
            return r

        def makeUnique(li):
            m={}
            for k in li:
                m[k]=0
            li=m.keys()
            li.sort()
            return li

        ise=makeUnique(isectIJ+isect)
        if not _pairsAreSorted(ise):
            print "ise",ise
        else:
            print ".",
        if len(ise)==0:
            try:
                del(self._TC[(i,j)])
                del(self._TC[(j,i)])
            except KeyError:
                pass
        else:
            self._TC[(i,j)]=ise
            self._TC[(j,i)]=ise
        #self._TC[(i,j)]=1
        #if len(isect)>0:
        #    print "isect",isect,isectIJ



    def extendTransitivePath2(self,i,j,k):
        """Extends the transitive path form i to j by k.

        This method is used in the inner most loop of the Floyd--Warshall
        algorithm in method TransitiveClosure(). The trick is to keep track
        of the exact mappings between aligned sites. Just tracking
        overlaps/doesn't overlap is not enough"""

        # REMOVE THIS:
        # needed for directed TC graph
        if i>j:
            try:
                self._TC[(i,j)]=self._TC[(j,i)]
            except KeyError:
                pass
        elif i==j:
            self._TC[(i,j)]=[99999]
            return

        if i==k or j==k:
            return 
        if j>k:
            self._TC[(j,k)]=self._TC[(k,j)]
        if i>k:
            try:
                self._TC[(i,k)]=self._TC[(k,i)]
            except KeyError:
                print self._TC.keys()
                raise
            

        try:
            isectIJ=self._TC[(i,j)]
        except KeyError:
            isectIJ=[]

        isectIK=self._TC[(i,k)]
        isectKJ=self._TC[(k,j)]

        OIS=isectIJ[:]

        # Find the intersecting sites.
        n,m=0,0
        while n<len(isectIK) and m<len(isectKJ):
            if isectIK[n][1]<isectKJ[m][0]:
                n+=1
            elif isectIK[n][1]>isectKJ[m][0]:
                m+=1
            else:
                assert(isectIK[n][1]==isectKJ[m][0])
                npair=(isectIK[n][0],isectKJ[m][1])
                try:
                    isectIJ.index(npair)
                except ValueError:
                    isectIJ.append((isectIK[n][0],isectKJ[m][1]))
                n,m=n+1,m+1

        isectIJ.sort()

        def _pairsAreSorted(l):
            r=1
            for i in range(1,len(l)):
                r=r and l[i-1][0]<l[i][0] and l[i-1][1]<l[i][1]
            return r
        try:
            assert(_pairsAreSorted(isectIJ))
        except AssertionError:
            pass
##            print "Inconsistent pair"
##            print "isectIJ (%d,%d)"%(i,j),OIS
##            print "isectIK (%d,%d)"%(i,k),isectIK
##            print "isectKJ (%d,%d)"%(k,j),isectKJ
##            print "new sorted isectIJ",isectIJ

##            print
##        if len(isectIJ)>0:
##            self._TC[(i,j)]=isectIJ
##            self._TC[(j,i)]=isectIJ
        self._TC[(i,j)]=isectIJ
        
                

    def TransitiveClosure(self):
        """Computes the transitive closure of the overlap graph.

        Uses the Floyd-Warshall algorithm for all-pairs-shortest-path
        problem with time complexity O(n^3). In theory, that can be
        improved to about O(n^2.3). The implementation is unnecessary
        strickt about having directed edges to both directions. This
        could, and should, be improved."""


        self._TC=self._overlaps.copy()  # Transitive closure
        n=len(self._pairs)
        if n<=2:
            return

        #Make non directed (i.e. edges both directions)
        for i,j in self._TC.keys():
            i,j=min((i,j)),max((i,j))
            self._TC[(j,i)]=self._TC[(i,j)]
        
        #print "BEFORE",self._TC.keys()
        def koe():
            for i,j in self._TC.keys():
                try:
                    assert(self._TC[(j,i)]==self._TC[(i,j)])
                except AssertionError:
                    print i,j,self._TC[(i,j)]
            
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    
                    if self._TC.has_key((i,k)) and self._TC.has_key((k,j)):
                        #print i,j,k
                        self.extendTransitivePath(i,j,k)
                        #self._TC[(i,j)]=1
                        #koe()
        #print "AFTER",self._TC.keys()
        koe()

##        for i in range(n):
##            l=len(self._pairs[i].sites)
##            self._TC[(i,i)]=zip(range(l),range(l))
##        print self._TC.keys()
        for p in self._TC.keys():
            print p,self._TC[p]
        #raise "EKA TEHTY"

    def overlapping(self,i,j):
        """Return true if i and j are overlapping (or same)"""
        return i==j or self._overlaps.has_key((i,j))


    def InconsistencyGraph2(self):
        """Computes the inconsistency graph for this multi module.

        The algorithm is more or less the same as described in
        Prohaska et.al."""

        self._inconsist={}

        def _DFSincons(i,trace):
            """Recursive part of the depth first search.

            Uses the function local variable dontFindOther which is a pair
            whose sites we should not find on our search.

            found is the list of found incompatible pairs."""
            if not self._overlapsM.has_key(i): return
            for j in self._overlapsM[i].keys():
                nPair=self._pairs[j]
                ntrace=nPair.makeTrace(trace)
                if not ntrace:
                    continue
                elif ntrace.isSubOrNoCommon2(dontFindOther):
                    _DFSincons(j,ntrace)
                else:
                    found.append(j)


        # Let's do the search starting from the end. This way we avoid
        # looping.
        for i in range(len(self._pairs)-1,0,-1):
            dontFindOther=self._pairs[i]
            found=[]
            _DFSincons(i,dontFindOther)
            for j in found:
                self._inconsist[(i,j)]=1
                self._inconsist[(j,i)]=1



    def InconsistencyGraph3(self):
        """Computes the inconsistency graph for this multi module.

        The algorithm is more or less the same as described in
        Prohaska et.al. This doesn't try fancy 'half matrix' things."""

        self._inconsist={}
        used={}
        def _DFSincons(i,trace):
            """Recursive part of the depth first search.

            Uses the function local variable dontFindOther which is a pair
            whose sites we should not find on our search.

            found is the list of found incompatible pairs."""
            goTo=[ x for (x,y) in self._overlaps.keys() if y==i]
            goTo=[ x for x in goTo if not used.has_key(x)]
            #print "goto",goTo

            used[i]=1
            for j in goTo:
                #if used.has_key(j): continue
                nPair=self._pairs[j]
                ntrace=nPair.makeTrace(trace)
                #print j,used.keys()
                #print "\n".join([x[:80] for x in str(ntrace).split("\n")])
                if ntrace:
                    if not ntrace.isSubOrNoCommon2(dontFindOther):
                        if not found.has_key(j):
                            found[j]=0
                        found[j]+=1
                    _DFSincons(j,ntrace)
                    #print "Found",j,"after",used.keys()
                    #pdb.set_trace()
            del(used[i])
            
        #print "\nStarting new IC graph"
        for i in range(len(self._pairs)):
            dontFindOther=self._pairs[i]
            found={}
            #print "Dont find other than:"
            #print "\n".join([x[:80] for x in str(dontFindOther).split("\n")])
            print "Looking for inconsistencies for",i
            superI=i
            _DFSincons(i,dontFindOther)
            print "found",found.keys(),used.keys()
            for j in found.keys():
                self._inconsist[(i,j)]=1
                self._inconsist[(j,i)]=1


    def InconsistencyGraph(self):
        """Computes the inconsistency graph for this multi module.

        This part too unnecessarily assumes directed graphs."""

        self._inconsist={}
        if len(self._pairs)<=2:
            return
        for i,j in self._TC.keys():
            isect=0
            for n in self._pairs[i].sNames:
                if self._pairs[j].DNApos.has_key(n):
                    isect+=1
            if isect==2 or (isect==1 and not self.overlapping(i,j)):
                self._inconsist[(i,j)]=1


    def connected(self,i,j):
        """Checks whether i and j are connected in the consistency graph."""
        return i==j or not self._inconsist.has_key((i,j))

    def cliqueToDot(self,j,strbuf=None):
        if not strbuf:
            strbuf=StringIO()
            
        clq=self._cliques[j]
        N=len(self._pairs)
        strbuf.write('subgraph cluster%d {\nlabel="Clique %s";\n%s;\n\n}'%(j,j,";\n".join(["n%dC%d"%(x,j) for x in clq])))
        for i in range(N):
            try:
                strbuf.write('n%dC%d [label="%s-%s CM%s"];\n'%(i,j,self._pairs[i].sNames[0][:6],self._pairs[i].sNames[1][:6],str(self._pairs[i].id[-3:])))
            except AttributeError:
                pass
            for k in range(i+1,N):
                if i!=k and self._overlaps.has_key((i,k)):
                    strbuf.write("n%dC%d -- n%dC%d;\n"%(i,j,k,j))
                if not self.connected(i,k):
                    strbuf.write("n%dC%d -- n%dC%d [color=red];\n"%(i,j,k,j))
        #strbuf.write("}\n")
        return strbuf

    
    def cliquesToDot(self):
        strbuf=StringIO()
        strbuf.write("graph G {\n")
        N=len(self._pairs)
        if hasattr(self,"_cliques"):
            for j in range(len(self._cliques)):
                self.cliqueToDot(j,strbuf)
        strbuf.write("}\n")
        return strbuf.getvalue()

    def inconsistToDot(self):
        """Returns the inconsistency graph in dot format.

        Dot is an graph layout program which is also known as graphviz.
        It is available for free from AT&T."""
        strbuf=StringIO()
        strbuf.write("graph G {\n")
        N=len(self._pairs)
        
        for i in range(N):
            try:
                strbuf.write('n%d [label="%s-%s CM%s"];\n'%(i,self._pairs[i].sNames[0][:3],self._pairs[i].sNames[1][:3],str(self._pairs[i].id[-3:])))
            except AttributeError:
                pass
            for j in range(i,N):
                if i!=j and self.connected(i,j):
                    strbuf.write("n%d -- n%d;\n"%(i,j))

        strbuf.write("}\n")
        return strbuf.getvalue()
                                 

    def allMaxCliques(self):
        """Compute all maximal cliques in the complement of the inconsistency graph.

        This uses the Bron, Kerbosch (CACM 1973: 16(9) p.575-7) algorithm
        to solve the NP complete problem (in worst case exponential time)"""

        N=len(self._pairs)

        if N<=2:
            self._cliques=[range(N)]
            return

        def extend2(old,ne,ce):
            """Extend procedure from Bron and Kerbosch"""
            #print old,ne,ce
            new=range(ce)
            minnod,nod=ce,0

            pos=ne # Default: Take the next candidate
            # Determine each counter value and look for minimum
            for i in range(ce):
                if minnod==0:
                    break
                p,count=old[i],0

                # Count disconnections
                for j in range(ne,ce):
                    if count >= minnod:
                        break
                    if not self.connected(p,old[j]):
                        count+=1
                        # Save position for potential candidate
                        pos=j

                # Test new minimum
                if count<minnod:
                    fixp=p
                    if i<ne:
                        try:
                            s=pos
                        except Exception:
                            print count,"<",minnod,old
                            raise
                    else:
                        s=i
                        nod=1
                    minnod=count
                    # Note: If fixed point initially chosen from candidates
                    # then number of disconnections will be preincreased by
                    # one

            # Backtrack cycle
            for nod in range(minnod+nod,0,-1):
                # Interchange
                old[s],old[ne]=old[ne],old[s]
                sel=old[ne]
                # Fill new set NOT
                newne=0
                for i in range(ne):  ## KJP: Make this with filter()
                    if self.connected(sel,old[i]):
                        new[newne]=old[i]
                        newne+=1
                        #print "n%d n%d connected, newne=%d"%(sel,old[i],newne)
                    #else:
                        #print "n%d n%d not connected ne=%d"%(sel,old[i],ne)
                # Fill new set CAND
                newce=newne
                for i in range(ne+1,ce):
                    if self.connected(sel,old[i]):
                        new[newce]=old[i]
                        newce=newce+1
                # Add to COMPSUB
                compsub[self.c]=sel;
                self.c+=1

                if newce==0:
                    self._cliques.append(compsub[:self.c])
                elif newne<newce:
                    try:
                        extend2(new,newne,newce)
                    except Exception:
                        print "new:",new,newne,newce
                        print "old:",old,ne,ce
                        raise
                # Remove from COMPSUB
                self.c-=1
                # Add to NOT
                ne+=1
                if nod>1:
                    # select a candidate disconnected to the fixed point
                    s=ne
                    # look for candidate
                    while self.connected(fixp,old[s]):
                        s+=1
                    try:
                        assert(ne<=s and not self.connected(fixp,old[s]))
                    except AssertionError:
                        print ne,"<=",s,self.connected(fixp,old[s])
                        raise
                        
        compsub=[-1]*N

        self._cliques=[]
        self.c=0
        extend2(range(N),0,N)


    def makeAllGreedyAligns(self,soCalledAll=20):
        """Runs the makeGreedyAlign2 for all seeds"""

        myPairs=self._pairs[:]
        myPairs.sort(lambda x,y:cmp(x.score,y.score))
        #a=[x.code for x in myPairs]
        #a.sort()
        print "AllGreedy",len(myPairs)
        multiMods=[]
        myPairsLim=max(len(myPairs)-soCalledAll,0)
        while len(myPairs)>myPairsLim:
            seed=myPairs.pop()
            if not seed.included:
                print "Seed %d:%g"%(seed.code,seed.score)
                multiMods.append(self.makeGreedyAlign2(seed))
                print "Result: %g"%(multiMods[-1].score)
        return multiMods        

    def makeGreedyAlign2(self,seed=None):
        """Makes a greedy local alignment from this module.

        Quite yacky method. Add overlapping pairwise alignments."""

        if not seed:
            # Decreasing order of pairs.
            myPairs=self._pairs[:]
            myPairs.sort(lambda x,y:cmp(x.score,y.score))
            seed=myPairs.pop() # The best scoring.
            
        mult={}
        sites=[]
        colCount=[]

        # Make the initial empty vectors
        for pair in self._pairs:
            for s in pair.sNames:
                mult[s]=[]


        alnLen=0
        Score=0.0
        pairCount=0

        class transact:
            def __init__(self):
                self.mult={}
            pass

        rollback=0
        
        pairMap={seed.code:seed}
        myPairs=pairMap.items()
        while len(myPairs)>0:
            #print "myPairs",len(myPairs)
            code,pair=myPairs.pop()
            del(pairMap[code])
            # Add pairs that intersects.



            seq1,seq2=pair.sNames

            # Prepare to back down from adding this pair
            back=transact()
            back.Sites=sites[:]
            back.ColCount=colCount[:]
            back.Score=Score
            back.pairCount=pairCount

            Score+=pair.score
            pairCount+=1
            for k,v in mult.items():
                back.mult[k]=v[:]


            for s in range(len(pair.sites)):
                p=0

                # Find new position
##                while p<alnLen and mult[seq1][p]<pair.sitePos[seq1][s] and \
##                      mult[seq2][p]<pair.sitePos[seq2][s]:
                while p<alnLen and (mult[seq1][p]<pair.sitePos[seq1][s] and \
                      mult[seq2][p]<pair.sitePos[seq2][s]):
                    p+=1

                # Verify location:
                p1,p2=p,p
                while p1<alnLen and mult[seq1][p1]==-1:
                    p1+=1
                while p2<alnLen and mult[seq2][p2]==-1:
                    p2+=1
                # Do not add this pair if the alignment would get busted.
                #if (p1<alnLen and (pair.sitePos[seq1][s]+pair.sites[s].width)>mult[seq1][p1]) or \
                #   (p2<alnLen and (pair.sitePos[seq2][s]+pair.sites[s].width)>mult[seq2][p2]):
                if (p1<alnLen and (pair.sitePos[seq1][s])>mult[seq1][p1]) or \
                   (p2<alnLen and (pair.sitePos[seq2][s])>mult[seq2][p2]):
                    mult=back.mult.copy()
                    colCount=back.ColCount
                    sites=back.Sites
                    Score=back.Score
                    pairCount=back.pairCount
                    alnLen=len(sites)
                    rollback=1
                    break
                
                
                # New column if needed
                if p==alnLen or ( (mult[seq1][p]>pair.sitePos[seq1][s] or \
                                   mult[seq1][p]==-1) and 
                                  (mult[seq2][p]>pair.sitePos[seq2][s] or \
                                   mult[seq2][p]==-1) ):
                    colCount.insert(p,0)
                    sites.insert(p,pair.sites[s])
                    alnLen+=1
                    for key in mult.keys():
                        mult[key].insert(p,-1)
                        try:
                            assert(alnLen==len(sites)==len(colCount)==len(mult[key]))
                        except AssertionError:
                            pdb.set_trace()
                try:
                    assert(mult[seq1][p] in (-1,pair.sitePos[seq1][s]))
                    assert(mult[seq2][p] in (-1,pair.sitePos[seq2][s]))
                except KeyError:
                    pdb.set_trace()
                    #print "Fitting %s (%g)"%(pair.id[-10:],pair.score)
                except AssertionError:
                    mult=back.mult.copy()
                    colCount=back.ColCount
                    sites=back.Sites
                    Score=back.Score
                    pairCount=back.pairCount
                    rollback=1
                    alnLen=len(sites)
                    #print "AE: non fitting alignment %s (%g)"%(pair.id[-10:],pair.score)
                    break

                # Set the values
                assert(mult[seq1][p] in (-1,pair.sitePos[seq1][s]))
                assert(mult[seq2][p] in (-1,pair.sitePos[seq2][s]))
                mult[seq1][p]=pair.sitePos[seq1][s]
                mult[seq2][p]=pair.sitePos[seq2][s]
                colCount[p]+=1

            pair.used=seed.code
            if rollback==0:
                pair.included=1
                [pairMap.setdefault(k,self._pairs[k]) for k in pair.overlaps if  abs(self._pairs[k].used)!=seed.code]
##            #myPairs.extend([self._pairs[x] for x in pair.overlaps])

##              Remove pairs that are already used.
##            #myPairs=[x for x in myPairs if not hasattr(x,"used")]

                myPairs=pairMap.items()
##            # Sort then so that myPairs.pop() is the best.
                myPairs.sort(lambda x,y:cmp(x[1].score,y[1].score))
            else:
                rollback=0

                


        return MultiOutModule(mult,sites,colCount,-1,0,Score,pairCount)

    def makeGreedyAlign(self):
        """Makes a greedy local alignment from this module.

        Very yacky method. Always add the best scoring, not yet added
        pairwise alignment."""

        # Decreasing order of pairs.
        myPairs=self._pairs[:]
        myPairs.sort(lambda x,y:cmp(y.score,x.score))
        mult={}
        back={}
        sites=[]
        colCount=[]

        # Make the initial empty vectors
        for pair in myPairs:
            for s in pair.sNames:
                mult[s]=[]


        alnLen=0
        Score=0.0
        for pair in myPairs:
            seq1,seq2=pair.sNames
            backSites=sites[:]
            backColCount=colCount[:]
            backScore=Score
            Score+=pair.score
            #Score+=1.0
            for k,v in mult.items():
                back[k]=v[:]
                try:
                    assert(alnLen==len(sites)==len(colCount)==len(mult[k]))
                except AssertionError:
                    print alnLen,len(sites),len(colCount),len(mult[k])
                    raise
            for s in range(len(pair.sites)):
                p=0

                # Find new position
##                while p<alnLen and mult[seq1][p]<pair.sitePos[seq1][s] and \
##                      mult[seq2][p]<pair.sitePos[seq2][s]:
                while p<alnLen and (mult[seq1][p]<pair.sitePos[seq1][s] and \
                      mult[seq2][p]<pair.sitePos[seq2][s]):
                    p+=1

                # Verify location:
                p1,p2=p,p
                while p1<alnLen and mult[seq1][p1]==-1:
                    p1+=1
                while p2<alnLen and mult[seq2][p2]==-1:
                    p2+=1
                # Do not add this if the alignment would get busted.
                #if (p1<alnLen and (pair.sitePos[seq1][s]+pair.sites[s].width)>mult[seq1][p1]) or \
                #   (p2<alnLen and (pair.sitePos[seq2][s]+pair.sites[s].width)>mult[seq2][p2]):
                if (p1<alnLen and (pair.sitePos[seq1][s])>mult[seq1][p1]) or \
                   (p2<alnLen and (pair.sitePos[seq2][s])>mult[seq2][p2]):
                    mult=back.copy()
                    colCount=backColCount
                    sites=backSites
                    Score=backScore
                    alnLen=len(sites)
                    break
                
                
                # New column if needed
                if p==alnLen or ( (mult[seq1][p]>pair.sitePos[seq1][s] or \
                                   mult[seq1][p]==-1) and 
                                  (mult[seq2][p]>pair.sitePos[seq2][s] or \
                                   mult[seq2][p]==-1) ):
                    colCount.insert(p,0)
                    sites.insert(p,pair.sites[s])
                    alnLen+=1
                    for key in mult.keys():
                        mult[key].insert(p,-1)
                        try:
                            assert(alnLen==len(sites)==len(colCount)==len(mult[key]))
                        except AssertionError:
                            pdb.set_trace()
                try:
                    assert(mult[seq1][p] in (-1,pair.sitePos[seq1][s]))
                    assert(mult[seq2][p] in (-1,pair.sitePos[seq2][s]))
                    #print "Fitting %s (%g)"%(pair.id[-10:],pair.score)
                except AssertionError:
                    mult=back.copy()
                    colCount=backColCount
                    sites=backSites
                    Score=backScore
                    alnLen=len(sites)
                    #print "AE: non fitting alignment %s (%g)"%(pair.id[-10:],pair.score)
                    break

                # Set the values
                assert(mult[seq1][p] in (-1,pair.sitePos[seq1][s]))
                assert(mult[seq2][p] in (-1,pair.sitePos[seq2][s]))
                mult[seq1][p]=pair.sitePos[seq1][s]
                mult[seq2][p]=pair.sitePos[seq2][s]
                colCount[p]+=1

            #print str(MultiOutModule(mult,sites,colCount,-1))

        return MultiOutModule(mult,sites,colCount,-1,0,Score)

    def makeMultiAlign(self,clqId):
        """Makes a dictionary of arrays representing the multiple alignment
        from clique clqId."""
        mult={}
        sites=[]
        colCount=[]
        clq=self._cliques[clqId]  # Operate on this clique
        clq.sort()
        # Make the initial empty vectors
        for i in clq:
            for s in self._pairs[i].sNames:
                mult[s]=[]

        N=len(clq)
        alnLen=0

        if clqId==0:
            print "#Clique 0 members"
            for c in clq:
                print str(self._pairs[c])
        
        for i in range(N):
            for j in range(i+1,N):
                # Overlaping pairs p1,p2
                p1,p2=self._pairs[clq[i]],self._pairs[clq[j]]
                try:
                    # Find out overlapping and non overlaping sequences
                    common=self._overlaps[(clq[j],clq[i])]
                    if p1.sNames[0]==common:
                        seq1=p1.sNames[1]
                    else:
                        seq1=p1.sNames[0]

                    if p2.sNames[0]==common:
                        seq2=p2.sNames[1]
                    else:
                        seq2=p2.sNames[0]
                        
                    # Overlapping sites s1 s2
                    p=0
                    try:
                        a,b=self._overlaps[(clq[i],clq[j])][0]
                    except ValueError:
                        print "Error coming"
                        print clq
                        print i,j
                        print self._overlaps
                        print self._overlaps[(clq[i],clq[j])]
                    except IndexError:
                        pass
                    goods,fails=0,0
                    for s1,s2 in self._overlaps[(clq[i],clq[j])]:
                        # Scan to correct position
                        try:
                            p<alnLen and mult[seq1][p]<p1.sitePos[seq1][s1] and \
                              mult[seq2][p]<p2.sitePos[seq2][s2] and \
                              mult[common][p]<p2.sitePos[common][s2]
                        except IndexError:
                            print mult
                            print p
                        while p<alnLen and mult[seq1][p]<p1.sitePos[seq1][s1] and \
                              mult[seq2][p]<p2.sitePos[seq2][s2] and \
                              mult[common][p]<p2.sitePos[common][s2]:
                            p+=1

                        try:
                            p==alnLen or ( (mult[seq1][p]>p1.sitePos[seq1][s1] or mult[seq1][p]==-1) and \
                              (mult[seq2][p]>p2.sitePos[seq2][s2] or mult[seq2][p]==-1) and \
                              (mult[common][p]>p2.sitePos[common][s2] or mult[common][p]==-1))
                        except IndexError:
                            traceback.print_exc()
                            print  "mult",mult
                            print "p=",p
                            print "p1",p1.sitePos
                            print "s1",s1
                            print "p2",p2.sitePos
                            print "s2",s2

                            
                        # Add new column if needed
                        if p==alnLen or ( (mult[seq1][p]>p1.sitePos[seq1][s1] or mult[seq1][p]==-1) and \
                              (mult[seq2][p]>p2.sitePos[seq2][s2] or mult[seq2][p]==-1) and \
                              (mult[common][p]>p2.sitePos[common][s2] or mult[common][p]==-1)):
                            for key in mult.keys():
                                mult[key].insert(p,-1)
                            colCount.insert(p,0)
                            sites.insert(p,p1.sites[s1])
                            alnLen+=1

                        #print seq1,s1,common,s2,seq2,s2,(clq[i],clq[j])
                        #print "mult[%s]=%s"%(seq1,str(mult[seq1]))
                        #print "mult[%s]=%s"%(common,str(mult[common]))
                        #print "mult[%s]=%s\n"%(seq2,str(mult[seq2]))
                        try:
                            assert(p1.sites[s1].match(p2.sites[s2]))
                            assert(mult[seq1][p]==-1 or mult[seq1][p]==p1.sitePos[seq1][s1])
                            assert(mult[common][p]==-1 or mult[common][p]==p1.sitePos[common][s1])
                            assert(mult[seq2][p]==-1 or mult[seq2][p]==p2.sitePos[seq2][s2])
                            goods+=1
                        except AssertionError:
                            fails+=1
                            #pdb.set_trace()
                            #continue
                            #return None
                            #raise
                            def pos_left(arr,x):
                                for p in range(len(arr)):
                                    if arr[p]>=0 and arr[p]>x:
                                        return p
                                return len(arr)
                            print "Graph size:",len(self._pairs)
                            print "CLIQUE:",clq,"Doing pairs",clq[i],clq[j]
                            print str(self._pairs[clq[i]])
                            print str(self._pairs[clq[j]])
                            #for c in clq:
                            #    print self._pairs[c].id
                            a=StringIO()
                            a.write("graph G {\n")
                            self.cliqueToDot(clqId,a)
                            a.write("}\n")
                            open("assertBreak%d.dot"%(clqId),"w").write(a.getvalue())
                            print "s1 %d,s2 %d, p %d"%(s1,s2,p)
                            print "p1.sites[%d]=%s"%(s1,str(p1.sites[s1]))
                            print "p2.sites[%d]=%s"%(s2,str(p2.sites[s2]))
                            print "mult[%s]=%s"%(seq1,str(mult[seq1]))
                            print "mult[%s]=%s"%(common,str(mult[common]))
                            print "mult[%s]=%s\n"%(seq2,str(mult[seq2]))
                            print "p1.sitePos[%s][%d]=%s"%(seq1,s1,str(p1.sitePos[seq1][s1]))
                            print "p1.sitePos[%s][%d]=%s"%(common,s1,str(p1.sitePos[common][s1]))
                            print "p2.sitePos[%s][%d]=%s"%(seq2,s2,str(p2.sitePos[seq2][s2]))
                            
                            mult[seq1].insert(pos_left(mult[seq1],p1.sitePos[seq1][s1]),"*")
                            mult[seq2].insert(pos_left(mult[seq2],p2.sitePos[seq2][s2]),"*")
                            mult[common].insert(pos_left(mult[common],p1.sitePos[common][s1]),"*")

                            for seq in mult.keys():
                                print "%10s %s"%(seq[:10],mult[seq])
                            #raise
                            #continue
                            return None
                        # (Re)set the values
                        mult[seq1][p]=p1.sitePos[seq1][s1]
                        mult[common][p]=p1.sitePos[common][s1]
                        mult[seq2][p]=p2.sitePos[seq2][s2]
                        colCount[p]+=1
                    if fails>0:
                        print "Good/fail= %d / %d"%(goods,fails)
                        print clq[i],clq[j]
                        print str(self._pairs[clq[i]])
                        print str(self._pairs[clq[j]])
                                  
                    #for key,val in mult.items():print "%s: %s"%(key[:5],str(val))
                    #print "      ",colCount
                        
                except KeyError:
                    pass
                    #print clq[i],clq[j]
                except ValueError:
                    print self._overlaps[(clq[i],clq[j])]
                    raise

        return MultiOutModule(mult,sites,colCount,clqId)

        
class Site:
    """Represents a binding site on the alignment"""
    def __init__(self,name,score,strand,width):
        self.name=name
        self.score=float(score)
        self.strand=strand
        self.width=int(width)

    def match(self,other):
        """Return true if two sites are for the same Transcription Factor"""
        return self.name==other.name and self.strand==other.strand and self.width==other.width

    def __str__(self):
        return "%s %s %d %g"%(self.name,self.strand,self.width,self.score)

    def toGFF(self,seq,method,spos):
        """Returns the site in GFF row format."""
        return "%s\t%s\t%s\t%s\t%d\t%d\t%g\t%s\t.\t"(seq,method,self.name,spos,spos+self.width,self.score,self.strand)
    
    def __repr__(self):
        return "Site(%s,%g,%s,%d)"%(self.name,self.score,self.strand,self.width)
class MultipleAlignment:
    """Main class for multiple alignment.

    Instance of this class is used to produce multiple alignment
    of binding site information. The input is any number of local
    pairwise alignments produced by mabs align."""
    
    def __init__(self):
        "Minimal constructor"
        self._modules=[]
        self._attrPat=re.compile(r"([A-Za-z][A-Za-z0-9_]+) (.*)")
        self._multiAlign=[]
        self.totLen=0

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
                if cisModRows<1:
                    if currMod:
                        self._modules.append(currMod)
                    lineModId="%s.%s.%d"%(fname,attribs["CM"],len(self._modules))
                    currMod=PairModule(lineModId,float(score),seq,int(start),int(stop),\
                                       seq2,start2,stop2)
                    cisModRows+=1
            else:
                cisModRows=0
                #assert(lineModId==currModId)
                siteId=Site(feat,float(score),strand,int(stop)-int(start))
                currMod.appendSite({seq:int(start)},siteId)
        if currMod:
            self._modules.append(currMod)

    def addGFFfile(self,fname):
        """Add a GFF file as input for the multiple alignment.

        The given GFF file should follow the format produced
        by mabs align command."""
        self._parseGFFfile(fname)
                
    def giveSequences(self):
        """Return information about sequences in the added GFF files.

        The return value is a list of tuples (seq,c) where seq is
        the name of the sequence and c is the number of cis modules
        found in the added GFF files for that sequence."""
        counts={}
        for module in self._modules:
            for seq in module.sNames:
                if not counts.has_key(seq):
                    counts[seq]=0
                counts[seq]=counts[seq]+1
        return counts.items()
    
    def removeSequences(self,sname):
        """Remove the sequence sname from the multiple alignment about
        to be produced."""
        self._modules=filter(lambda x:not x.sNames.has_key(sname),self._modules)


    def __len__(self):
        return len(self.resAligns)

    def __getitem__(self,i):
        try:
            return self.resAligns[i]
        except AttributeError:
            raise IndexError("Index out of bounds")

    
    def multiAlign(self):
        """Compute the local multiple alignment.

        The input and output is done through instance of this class with
        other methods."""

        self.resAligns=[]
        nAligns=[]
        import random
        #random.shuffle(self._modules)
        while len(self._modules)>0:
            nAligns=[]
            r=None
            p=self._modules.pop()
            for mmod in self._multiAlign:
                if mmod.addPair(p):
                    if not r:
                        r=mmod
                        nAligns.append(r)
                    else:
                        #print "Uniting two multi modules"
                        r.uniteButLast(mmod)
                else:
                    nAligns.append(mmod)
            if not r:
                nMulti=MultiModule()
                nMulti.addPair(p)
                nAligns.append(nMulti)

            self._multiAlign=nAligns[:]

##        print len(self._multiAlign),
        self._multiAlign=[x for x in self._multiAlign if len(x._pairs)>1 ]
##        print len(self._multiAlign)
##        print [len(x._pairs) for x in self._multiAlign]
        #self._multiAlign=filter(lambda x:len(x._pairs)>1,self._multiAlign)
        for mmod in self._multiAlign:
##            self.resAligns.append(mmod.makeGreedyAlign2())
            self.resAligns.extend(mmod.makeAllGreedyAligns())

        self.resAligns.sort(lambda x,y:cmp(y.score,x.score))

    def DOESNTWORK(self):
        assert(None)
        # Multiple alignment by Prohaska et.al.
        c=0
        for mmod in self._multiAlign:
            #mmod.TransitiveClosure()
            #mmod.InconsistencyGraph()
            mmod.InconsistencyGraph3()
            mmod.allMaxCliques()

            for i in range(len(mmod._cliques)):
                try:
                    #print "Module",c,"Cliques:",len(mmod._cliques),[x.id[-10:] for x in mmod._pairs]
                    self.resAligns.append(mmod.makeMultiAlign(i))
                except AssertionError:
                    print "AE"
                    print mmod._overlaps.keys()
                    raise
            c=c+1
            
    def __str__(self):
        strbuf=StringIO()
        for mod in self._modules:
            strbuf.write(str(mod))
        return strbuf.getvalue()
        

def apu(datat=1):
    a=MultipleAlignment()
    if datat==0:
        print "MYF5"
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hmmyf5.align.opt.gff")
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hrmyf5.align.opt.gff")
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/mrmyf5.align.opt.gff")
    elif datat==1:
        print "ENSall"
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/ENS.10ekaa.gff")
    elif datat==2:
        print "Multikoe"
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/mabslib/tmp/cliqueMin.gff")
    elif datat==3:
        print "bestmyf5"
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/mabslib/tmp/bestmyf5.gff")
    a.multiAlign()


    for x in a:
        print str(x)
##    for b in range(len(a._multiAlign)):
##        for c in range(len(a._multiAlign[b]._cliques)):
##            if len(a._multiAlign[b]._cliques[c])<2: continue
##            #print "Module",b,"clique",c
##            #print
##            try:
##                x=str(a._multiAlign[b].makeMultiAlign(c))
##                print x
##            except AssertionError:
##                raise
##                pass
##        open("tmp.%d.dot"%(b),"w").write(a._multiAlign[b].inconsistToDot())
##        #print "Module",b,"cliques",len(a._multiAlign[b]._cliques)

##    if datat==3:
##        goodOverlap=[(3, 2), (3, 1), (0, 1), (5, 1), (2, 7), (8, 4), (4, 1), (6, 1), (1, 6), (7, 2), (7, 1), (2, 1), (4, 8), (1, 0), (1, 3), (1, 2), (1, 5), (1, 4), (1, 7), (2, 3)]
##        goodOverlap.sort()
##        goodMap=['f5.gff.3.4', 'f5.gff.2.1', 'f5.gff.2.5', 'f5.gff.7.8', '5.gff.10.2', '5.gff.11.0', 'f5.gff.5.3', 'f5.gff.3.6', 'f5.gff.4.7']

##        b=a._multiAlign[0]
##        helpOverlap=[ (goodMap.index(b._pairs[x].id[-10:]),goodMap.index(b._pairs[y].id[-10:])) for (x,y) in b._overlaps.keys()]
##        helpOverlap.sort()
##        print helpOverlap
##        print goodOverlap
##        assert(helpOverlap==goodOverlap)

##        goodIncons=[(1, 5), (2, 4), (3, 4), (4, 2), (4, 3), (5, 1)]
##        helpIncons=[ (goodMap.index(b._pairs[x].id[-10:]),goodMap.index(b._pairs[y].id[-10:])) for (x,y) in b._inconsist.keys()]
##        helpIncons2=[ (b._pairs[x].id[-10:],b._pairs[y].id[-10:]) for (x,y) in b._inconsist.keys()]
##        helpIncons.sort()
##        print helpIncons2
##        print 
##        print helpIncons
##        print goodIncons
##        assert(helpIncons==goodIncons)

##    if datat==0:
##        print "Should be\n(cliques,inconsists,pairs,overlapsM) [(1, 0, 4, 4), (36, 94, 18, 50), (1, 0, 3, 3), (1, 0, 2, 1), (1, 0, 3, 2)]"
##    x=[ (len(x._cliques),len(x._inconsist),len(x._pairs),reduce(lambda z,y:z+y,map(len,x._overlapsM.values()))) for x in a._multiAlign]
##    #x.sort()
##    print "(cliques,inconsists,pairs,overlapsM)",x
####    for b in range(len(a)):
####        x=str(a[b])
####        print "Align %d\n"%(b),x
    

if __name__=="__main__":
    import profile,pstats,sys
    if 0:
        profile.run("apu(0)","mabs.prof")
        p=pstats.Stats("mabs.prof")
        p.sort_stats('time').print_stats(20)
    else:
        if len(sys.argv)>1:
            code=int(sys.argv[1])
        else:
            code=0
        pdb.run("apu(code)")
        #apu(code)
        
