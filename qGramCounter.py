from Bio.SeqIO import FASTA
import sys

file=FASTA.FastaReader(open(sys.argv[1]))


for rec in file:
    q=4
    grams={}
    print rec.id,len(rec.seq)
    for i in range(1,len(rec.seq)):
        if (i%1000)==0:
            print i,
            sys.stdout.flush()
        for j in range(1,min(i,q)+1):
            gram=rec.seq[i-j:i]
            if grams.has_key(gram):
                grams[gram]=grams[gram]+1
            else:
                grams[gram]=1
                
    print rec.id
