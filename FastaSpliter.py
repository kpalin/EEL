from Bio.SeqIO import FASTA

file=FASTA.FastaReader(open("Nmyc.Enhancers.50kBp.fa"))


for rec in file:
    fname="Nmyc.Enhancer.50kBp.%s.fa"%(rec.id)
    fout=FASTA.FastaWriter(open(fname,"w"))
    fout.write(rec)
    print rec.id
