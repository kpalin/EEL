#!/usr/bin/python2.2

from random import shuffle

a,c,g,t=0.2,0.3,0.3,0.2
n=50000

s=["A"]*int(n*a)+["C"]*int(n*c)+["G"]*int(n*g)+["T"]*int(n*t)
shuffle(s)
s="".join(s)

fout=open("random.fa","w")
fout.write(">Random Sequence with %.2f %.2f %.2f %.2f of A C G T\n"%(a,c,g,t))

lineLen=70

for i in range(0,len(s),lineLen):
    fout.write(s[i:i+lineLen]+"\n")

fout.close()
