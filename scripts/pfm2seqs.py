#!/usr/bin/python
# -*- coding: UTF-8 -*- 


#
# $Log$
#

# For converting pfm files to format accepted by ModuleFinder by Bulyk et.al.


import sys

from eellib.Matrix import Matrix

limit=float(sys.argv[1])
for i in range(2,len(sys.argv)):
    fname=sys.argv[i]
    try:
        m=Matrix(fname)
    except ValueError,e:
        sys.stderr.write("Couldn't read %s (%s)\n"%(fname,str(e)))
        pass
    if len(m)>15:
        sys.stderr.write("N%d %s\n"%(len(m),fname))
        sbl=m.seqsBetterThan(limit)
        #print "\n".join([x[0] for x in sbl])
