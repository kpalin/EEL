import sys

#
# $Log$
#

from Bio.SeqIO import FASTA

file=FASTA.FastaReader(open(sys.argv[1]))


count=0
for rec in file:
    fname="%s.%d.fa"%(sys.argv[1],count)
    count+=1
    fout=FASTA.FastaWriter(open(fname,"w"))
    fout.write(rec)
    print rec.id
