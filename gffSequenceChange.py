"""Convert GFF sequence names and coordinates.

usage:
python2.2 gffSequenceChange.py [-S sname] [-s tname] [-p offset] -f fname [-a]

-S sname    Change sequences with sname
-s tname    Change the sequence name to tname
-p offset   Add offset to all coordinates
-f fname    In file 'fname'
-a          Output all lines, not just changed ones.
"""
import sys,re,getopt


comopt,args = getopt.getopt(sys.argv[1:],"S:s:p:f:a")

FromChg,toChg,coordChg,fname,outputAll="","",0,"need.a.file.gff",0
for opt,val in comopt:
    if opt=="-S":
        fromChg=val
    elif opt=="-s":
        toChg=val
    elif opt=="-p":
        coordChg=int(val)
    elif opt=="-f":
        fname=val
    elif opt=="-a":
        outputAll=1



for line in map(lambda x:x.split("\t"),open(fname).readlines()):
    if line[0]==fromChg:
        line[0]=toChg
        line[3]=str(int(line[3])+coordChg)
        line[4]=str(int(line[4])+coordChg)
        print "\t".join(line),
    elif outputAll:
        print "\t".join(line),

