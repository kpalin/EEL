"""Python interface and logic for multiple alignment."""

#
# $Log$
#



class PairModule:
    """A pairwise cis module from a GFF file.

    Use as a representation of the pairwise alignment cis
    module from the mabs alignment GFF files. This class
    has several important functions for adding pairs to form
    multiple alignments."""

    def __init__(self,cisCode,sname1,spos1,epos1,sname2,spos2,epos2):
        """Constructor defining the sequence name and position
        of the cis module.

        cisCode is the unique code for this cis module.
        The parameter sname1 is the name for the first sequence,
        spos1 is the start and epos is the end position on the
        first sequence DNA. The parameters sname2, spos2 and epos2
        are the same for the second sequence."""
        self.id=cisCode

        self.sNames=[sname1,sname2]
        self.sNames.sort()
        
        self.DNApos={}
        self.DNApos[sname1]=(spos1,epos1)
        self.DNApos[sname2]=(spos2,epos2)

    def __cmp__(self,other):
        """Comparison between two PairModule:s

        The PairModules are ordered in lexicographic order with tuple
        spos1,epos1,spos2,epos2  where the sequences are in alphabetic
        order."""
        r=cmp(self.DNApos[self.sNames[0]],other.DNApos[other.sNames[0]])
        if r==0:
            r=cmp(self.DNApos[self.sNames[1]],other.DNApos[other.sNames[1]])
        if self.id==other.id:
            assert(self.sNames==other.sNames)
        return r

    def intersects(self,other):
        """Checks if there is intersection between the two PairModules.

        Return true if the modules share a common TF binding site,
        false otherwise."""
        pass

    def appendSite(self,sPos,site):
        """Append a new pair of sites belonging to this cis module.

        The parameter sPos is a map from sequence name to tuple of sites
        DNA coordinate and match score. site is the TF binding site
        identification."""
        pass



class MultiModule:
    """Class that forms to be the cis module in multiple sequences.

    Input is a long list of pair modules."""

    def __init__(self):
        """Minimal constructor.

        Most work is done while adding the pairwise alignments."""
        pass

    def addPair(self,pair):
        """Adds a PairModule to this multiple alignment.

        The multiple alignment is formed by aligning all pairs
        of the sequences in the alignment with each other. These
        pairwise alignments are done before the multiple alignment."""
        pass

    def enlargeAlignment(self,pair):
        """Adds a new sequence to the multiple alignment.

        Adds a new PairModule to the MultiModule such that
        the PairModule contains a sequence not yet in MultiModule.
        This gives a new sequence to the multiple alignment."""
        pass
    
    def strengthenAlignment(self,pair):
        """Adds a new alignment between pairs that are already
        in the multiple alignment.

        With this method we should add alignments between all
        sequences in the alignment."""
        pass
    

        
    

class MultipleAlignment:
    """Main class for multiple alignment.

    Instance of this class is used to produce multiple alignment
    of binding site information. The input is any number of local
    pairwise alignments produced by mabs align."""
    
    def __init__(self):
        "Minimal constructor"
        pass

    def addGFFfile(self,fname):
        """Add a GFF file as input for the multiple alignment.

        The given GFF file should follow the format produced
        by mabs align command."""
        pass

    def giveSequences(self):
        """Return information about sequences in the added GFF files.

        The return value is a list of tuples (seq,c) where seq is
        the name of the sequence and c is the number of cis modules
        found in the added GFF files for that sequence."""
        pass
    
    def removeSequences(self,sname):
        """Remove the sequence sname from the multiple alignment about
        to be produced."""
        pass


    def multiAlign(self):
        """Compute the local multiple alignment.

        The input and output is done through instance of this class with
        other methods."""
