
import sys,re,getopt


comopt,args = getopt.getopt(sys.argv[1:],"S:s:p:f:a")

outputAll=0
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

