"""Python interface and logic for multiple alignment."""

#
# $Log$
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
        for common in self.sNames:
            if other.DNApos.has_key(common):
                break
        try:
            assert(self.DNApos[common][0]<other.DNApos[common][0]<self.DNApos[common][1])
        except AssertionError,KeyError:
            return 0
        for (pos,sites) in zip(self.sitePos,self.sites):
            for p,s in zip(other.sitePos,other.sites):
                if p==pos and sites.match(s):
                    return 1
        return 0
                

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
                         


class MultiModule:
    """Class that forms to be the cis module in multiple sequences.

    Input is a long list of pair modules."""

    def __init__(self):
        """Minimal constructor.

        Most work is done while adding the pairwise alignments."""
        self._TFseqs={}
        self.sites=[]
        self._pairs={}
        pass


    def __str__(self):
        return str(self._pairs.items())

    def addPair(self,pair):
        """Adds a PairModule to this multiple alignment.

        The multiple alignment is formed by aligning all pairs
        of the sequences in the alignment with each other. These
        pairwise alignments are done before the multiple alignment."""
        r=0
        if len(self._pairs)==0:
            self._TFseqs=pair.sitePos.copy()
            self._sites=pair.sites[:]
            r=1
            for name in pair.sNames:
                if not self._pairs.has_key(name):
                    self._pairs[name]=[]
                self._pairs[name].append(pair)
        elif self._TFseqs.has_key(pair.sNames[0]) and self._TFseqs.has_key(pair.sNames[1]):
            r=self.strengthenAlignment(pair)
        elif self._TFseqs.has_key(pair.sNames[0]) or self._TFseqs.has_key(pair.sNames[1]):
            r=self.enlargeAlignment(pair)
        return r

    def enlargeAlignment(self,pair):
        """Adds a new sequence to the multiple alignment.

        Adds a new PairModule to the MultiModule such that
        the PairModule contains a sequence not yet in MultiModule.
        This gives a new sequence to the multiple alignment."""

        # Do we have a common sequence?
        haveKey=0
        for common in pair.sNames:
            if self._pairs.has_key(common):
                haveKey=1
                break
        if haveKey==0:
            return None

        # Do we have an overlapping pairwise alignment?
        intersects=0
        for p in self._pairs[common]:
            intersects=p.intersects(pair)
            if intersects:
                break
        if intersects==0:
            return None

        for name in pair.sNames:
            if not self._pairs.has_key(name):
                self._pairs[name]=[]
            self._pairs[name].append(pair)

        i=0
        TFseqs={}
        for k in self._TFseqs.keys():
            TFseqs[k]=0
        while i<len(self._sites):
            for k in TFseqs.keys():
                
            
        return self
    
    def strengthenAlignment(self,pair):
        """Adds a new alignment between pairs that are already
        in the multiple alignment.

        With this method we should add alignments between all
        sequences in the alignment."""
        seq1,seq2=pair.sNames
        for p in self._pairs[seq1]:
            intersects=p.intersects(pair)
            if intersects:
                break
        if intersects==0:
            return None
        intersects=0
        for q in self._pairs[seq2]:
            intersects=q.intersects(pair)
            if intersects:
                break
        if intersects==0:
            return None

        for name in pair.sNames:
            if not self._pairs.has_key(name):
                self._pairs[name]=[]
            self._pairs[name].append(pair)
        return self
    

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
        pass


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
            if len(line)<1 or line[0]=='#':  ## Skip empty and comment lines.
                continue
            ## Parse a GFF line
            line=line.strip()
            seq,src,feat,start,stop,score,strand,frame,attribs=line.split('#')[0].split("\t",8)
            attribs=self._parseAttribs(attribs)

            lineModId=fname+attribs["CM"]
            # Start of a new cis module
            if feat=='CisModule':
                seq2=attribs["Target"].strip('"')
                start2=int(attribs["Start"])
                stop2=int(attribs["End"])
                if currModId!=lineModId:
                    if currMod:
                        self._modules.append(currMod)
                    currMod=PairModule(lineModId,float(score),seq,int(start),int(stop),\
                                       seq2,start2,stop2)
                    currModId=lineModId
            else:
                assert(lineModId==currModId)
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


    def multiAlign(self):
        """Compute the local multiple alignment.

        The input and output is done through instance of this class with
        other methods."""
        
        # Add the pair modules in decreasing order of alignment score
        self._modules.sort(lambda x,y:cmp(x.score,y.score))

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

    def __str__(self):
        strbuf=StringIO()
        for mod in self._modules:
            strbuf.write(str(mod))
        return strbuf.getvalue()
        


if __name__=="__main__":
    a=MultipleAlignment()
    #a.addGFFfile("/home/kpalin/tyot/comparative/mabs/fly.abs9.align.gff")
    a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hmmyf5.align.opt.gff")
    a.addGFFfile("/home/kpalin/tyot/comparative/mabs/hrmyf5.align.opt.gff")
    #a.addGFFfile("/home/kpalin/tyot/comparative/mabs/mrmyf5.align.opt.gff")
    a._modules.sort()
    print str(a)
    print a.giveSequences()
    a.multiAlign()
