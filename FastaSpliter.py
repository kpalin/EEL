import sys

#
# $Log$
# Revision 1.2  2004/02/23 12:22:57  kpalin
# Updates for per gene orthologous runs.
#
#

from Bio.SeqIO import FASTA

file=FASTA.FastaReader(open(sys.argv[1]))

idName=None
if len(sys.argv)>2:
    idName=sys.argv[2]

count=0
for rec in file:
    if not idName:
        idName=str(count)
    fname="%s.%s.fa"%(rec.id,idName)
    count+=1
    fout=FASTA.FastaWriter(open(fname,"w"))
    fout.write(rec)
    print rec.id
