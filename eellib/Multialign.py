"""Python interface and logic for multiple alignment."""

#
# $Log$
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


import xreadlines,re
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
        
        self.DNApos={}
        self.DNApos[sname1]=(spos1,epos1)
        self.DNApos[sname2]=(spos2,epos2)
        self.sitePos={sname1:[],sname2:[]}
        self.sites=[]

    def __cmp__(self,other):
        """Comparison between two PairModule:s

        The PairModules are partially ordered along their alphabetically
        first common sequence."""
        haveCommon=1
        for common in self.sNames:
            if common in other.sNames:
                haveCommon=1
                break
        try:
            assert(haveCommon==1)
        except AssertionError:
            raise NotComparable(self,other)
            
        r=cmp(self.DNApos[common],other.DNApos[common])
        if r==0:
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
        for seq in self.sNames:
            self.sitePos[seq].sort()


    def __str__(self):
        self._finalize()
        
        strbuf=StringIO()
        
        cmFormat='%s\tmalign\tCisModule\t%d\t%d\t%g\t.\t.\tTarget "%s";\tStart %d;\tEnd %d;\tStrand .;\tFrame .;\tCM %s;\n'

        seq1,seq2=self.sNames
        strbuf.write(cmFormat%(seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.score,seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.id))
        strbuf.write(cmFormat%(seq2,self.DNApos[seq2][0],self.DNApos[seq2][1],self.score,seq1,self.DNApos[seq1][0],self.DNApos[seq1][1],self.id))

        return strbuf.getvalue()

    def toGFF(self):
        """Returns GFF representation of the pairModule."""
        self._finalize()
        
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
    def __init__(self,aligns,sites,colCounts,id):
        self.align=aligns
        self.sites=sites
        self.colCount=colCounts
        self.id=id
        self.n=len(self.sites)

        self.sb=StringIO()
        self.seqs=self.align.keys()
        self.seqContext=10
        self.writeSeqNames()
        self.writePosAln()
       
    def __len__(self):
        return len(self.sites)

    def writeSeqNames(self):
        self.sb.write("Sequences:\n")
        for seq in self.seqs:
            self.sb.write("%s\n"%(seq))
        self.sb.write("\n")

    def writePosAln(self):
        self.sb.write(" "*26+"  ".join(["%-10s"%(x[:10]) for x in self.seqs]))
        for i in range(self.n):
            self.sb.write("\n%-20s %s %2d "%(self.sites[i].name[-20:],self.sites[i].strand,self.colCount[i]))
            self.sb.write("  ".join(["%10d"%(self.align[x][i]) for x in self.seqs]))
        self.sb.write("\n")

    def __str__(self):
        return self.sb.getvalue()


    def writeAlnSeq(self,linelen=60):
        pos={}
        aln={}
        self.sb.write("\n")
        for (seq,i) in zip(self.seqs,range(len(self.seqs))):
            pos[seq]=self.align[seq][0]-self.seqContext
            aln[seq]=self.seqAln[seq].getvalue()
            self.sb.write("Sequence %d : %s\n"%(i,seq))

        alnLen=len(aln.values()[0])
        self.sb.write("\n")
        for p in range(0,alnLen,linelen):
            for seq in self.seqs:
                try:
                    alnLine=aln[seq][p:p+linelen]
                except KeyError:
                    print seq,self.seqs,aln.keys()
                    print self.sb.getvalue()
                    raise
                self.sb.write("%7d : %s\n"%(pos[seq],alnLine))
                pos[seq]=pos[seq]+linelen-alnLine.count("-")
            self.sb.write("\n")
            
                
            

    def strAln(self,fullSeqs):
        prevPos={}
        self.seqAln={}
        for seq in self.seqs:
            prevPos[seq]=self.align[seq][0]-self.seqContext
            self.seqAln[seq]=StringIO()

        #Give 10 bp of context before alignments.
        intervals=[self.seqContext]
        for p in range(1,self.n):
            intervals.append(0)
            for seq in self.seqs:
                try:
                    thisInt=self.align[seq][p]-(self.sites[p-1].width+self.align[seq][p-1])
                    intervals[-1]=max(intervals[-1],thisInt)
                except TypeError:
                    pass

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
        pass


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
            if not intersects:
                return 0

        for name in pair.sNames:
            if not self._pairsBySeq.has_key(name):
                self._pairsBySeq[name]=[]
            self._pairsBySeq[name].append((pair,n))
        self._pairs.append(pair)
        
        return 1                
        


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

        for k in range(n):
            for i in range(n):
                for j in range(n):
                    
                    if self._TC.has_key((i,k)) and self._TC.has_key((k,j)):
                        self._TC[(i,j)]=1

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
            if isect==2 or (isect==1 and not self._overlaps.has_key((i,j))):
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
            for k in range(i,N):
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
                    for s1,s2 in self._overlaps[(clq[i],clq[j])]:
                        # Scan to correct position
                        while p<alnLen and mult[seq1][p]<p1.sitePos[seq1][s1] and \
                              mult[seq2][p]<p2.sitePos[seq2][s2] and \
                              mult[common][p]<p2.sitePos[common][s2]:
                            p+=1

                        # Add new column if needed
                        if p==alnLen or ( (mult[seq1][p]>p1.sitePos[seq1][s1] or mult[seq1][p]==-1) and \
                              (mult[seq2][p]>p2.sitePos[seq2][s2] or mult[seq2][p]==-1) and \
                              (mult[common][p]>p2.sitePos[common][s2] or mult[common][p]==-1)):
                            for key in mult.keys():
                                mult[key].insert(p,-1)
                            colCount.insert(p,0)
                            sites.insert(p,p1.sites[s1])
                            alnLen+=1

                        try:
                            assert(p1.sites[s1].match(p2.sites[s2]))
                            assert(mult[seq1][p]==-1 or mult[seq1][p]==p1.sitePos[seq1][s1])
                            assert(mult[common][p]==-1 or mult[common][p]==p1.sitePos[common][s1])
                            assert(mult[seq2][p]==-1 or mult[seq2][p]==p2.sitePos[seq2][s2])
                        except AssertionError:
                            raise
                            def pos_left(arr,x):
                                for p in range(len(arr)):
                                    if arr[p]>=0 and arr[p]>x:
                                        return p
                                return len(arr)
                            print "Graph size:",len(self._pairs)
                            print "CLIQUE:",clq
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
                            raise
                        # (Re)set the values
                        mult[seq1][p]=p1.sitePos[seq1][s1]
                        mult[common][p]=p1.sitePos[common][s1]
                        mult[seq2][p]=p2.sitePos[seq2][s2]
                        colCount[p]+=1
                        
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
            lineModId=fname+attribs["CM"]
            # Start of a new cis module
            if feat=='CisModule':
                try:
                    seq2=attribs["Target"].strip('"')
                except TypeError:
                    seq2=attribs["Target"].replace('"','')
                start2=int(attribs["Start"])
                stop2=int(attribs["End"])
                if currModId!=lineModId:
                    if currMod:
                        self._modules.append(currMod)
                    currMod=PairModule(lineModId,float(score),seq,int(start),int(stop),\
                                       seq2,start2,stop2)
                    currModId=lineModId
            else:
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
        return self.resAligns[i]

    
    def multiAlign(self):
        """Compute the local multiple alignment.

        The input and output is done through instance of this class with
        other methods."""

        self.resAligns=[]
        while len(self._modules)>0:
            r=0
            p=self._modules.pop()
            for mmod in self._multiAlign:
                if mmod.addPair(p):
                    r=1
            if r==0:
                nMulti=MultiModule()
                nMulti.addPair(p)
                self._multiAlign.append(nMulti)
        # Multiple alignment by Prohaska et.al.
        for mmod in self._multiAlign:
            mmod.TransitiveClosure()
            mmod.InconsistencyGraph()
            mmod.allMaxCliques()

            for i in range(len(mmod._cliques)):
                if len(mmod._cliques[i])<2:continue
                try:
                    self.resAligns.append(mmod.makeMultiAlign(i))
                except AssertionError:
                    raise
            
    def __str__(self):
        strbuf=StringIO()
        for mod in self._modules:
            strbuf.write(str(mod))
        return strbuf.getvalue()
        

def apu(datat=1):
    a=MultipleAlignment()
    if datat==0:
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hmmyf5.align.opt.gff")
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hrmyf5.align.opt.gff")
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/mrmyf5.align.opt.gff")
    elif datat==1:
        a.addGFFfile("/home/kpalin/tyot/comparative/mabs/ENS.10ekaa.gff")
    a.multiAlign()

    return 
    for b in range(len(a._multiAlign)):
        for c in range(len(a._multiAlign[b]._cliques)):
            print "Module",b,"clique",c
            #print
            try:
                x=str(a._multiAlign[b].makeMultiAlign(c))
                print x
            except AssertionError:
                pass
        open("tmp.%d.dot"%(b),"w").write(a._multiAlign[b].cliquesToDot())
        print "Module",b,"cliques",len(a._multiAlign[b]._cliques)

    x=[len(x._cliques) for x in a._multiAlign]
    x.sort()
    print x
##    for b in range(len(a)):
##        x=str(a[b])
##        print "Align %d\n"%(b),x
    

if __name__=="__main__":
    import profile,pstats
    profile.run("apu(1)","mabs.prof")
    p=pstats.Stats("mabs.prof")
    p.sort_stats('time').print_stats(20)
