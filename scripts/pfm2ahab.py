#!/usr/bin/python
# -*- coding: UTF-8 -*- 


#
# $Log$
#

# For converting pfm files to format accepted by ahab (by rajewsky vergassola gaul siggia)


import sys

from eellib.Matrix import Matrix


for i in range(1,len(sys.argv)):
    fname=sys.argv[i]
    sys.stderr.write("N%d %s\n"%(i,fname))
    print Matrix(fname).toAhab(i)
