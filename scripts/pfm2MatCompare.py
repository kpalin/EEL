#!/usr/bin/python
# -*- coding: UTF-8 -*- 


#
# $Log$
# Revision 1.1  2005/01/25 09:11:49  kpalin
# An assist script to convert pfm files to format uploadable to ahab.
#
#

# For converting pfm files to format accepted by ahab (by rajewsky vergassola gaul siggia)


import sys

from eellib.Matrix import Matrix
from StringIO import StringIO

for i in range(1,len(sys.argv)):
    fname=sys.argv[i]
    sys.stderr.write("M%d %s\n"%(i,fname))
    try:
        print Matrix(fname).toMatCompare().getvalue()
    except ValueError,e:
        sys.stderr.write("Couldn't read %s (%s)\n"%(fname,str(e)))
        pass
