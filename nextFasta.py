from xmlrpclib import ServerProxy


#
# $Log$
#

import sys

servName="alkokrunni.cs.helsinki.fi"
servPort=80


server=ServerProxy("http://%s/~kpalin/mabsCommandServerCGI.py"%(servName))


import sys
import socket

try:
    cID="%s:%s"%(socket.gethostname(),sys.argv[1])
    print server.nextFasta(cID)
except Exception,e:
    print >> sys.stderr, e
    print "xTryAgain"
