"""Converts a GFF file to chromosomal coordinates.

GFF file with sequence names like HUMAN|1.12345-19999
are converted to names like 1 with the start coordinate
added to the feature coordinates. Just remember to properly GREP
the input before running this!!.

Output comes to standard output. Filename is given on commandline."""
import sys,re

##
## $Log$
## 


fname=sys.argv[1]
import re

coordReg=re.compile(r'\w+\|(\d+).(\d+)-(\d+)')

for line in map(lambda x:x.split("\t",8),open(fname).readlines()):
    coordMatch=coordReg.search(line[0])
    if coordMatch:
        start,end=coordMatch.span()
        line[0]="%s%d%s"%(line[0][:start],int(coordMatch.group(1)),line[0][end:])
        coordChg=int(coordMatch.group(2))
        assert(coordChg<int(coordMatch.group(3)))
        line[3]=str(int(line[3])+coordChg)
        line[4]=str(int(line[4])+coordChg)
    coordMatch=coordReg.search(line[-1])
    if coordMatch:
        start,end=coordMatch.span()
        line[-1]="%s%d%s"%(line[-1][:start],int(coordMatch.group(1)),line[-1][end:])
        start=line[-1].index("Start ")+6
        end=line[-1].index(";",start)
        line[-1]="%s%d%s"%(line[-1][:start],int(coordMatch.group(2))+int(line[-1][start:end]),line[-1][end:])
        start=line[-1].index("End ")+4
        end=line[-1].index(";",start)
        line[-1]="%s%d%s"%(line[-1][:start],int(coordMatch.group(2))+int(line[-1][start:end]),line[-1][end:])
    print "\t".join(line),

