#!/usr/bin/python2.2

import sys

line=sys.stdin.readline()
while len(line)>0:
    if line[0]=="#":
        print line,
    else:
        parts=list(line.split("\t"))
        parts[6]="."
        print "\t".join(parts),
    line=sys.stdin.readline()



