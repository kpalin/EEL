#!/usr/bin/python2.2
import socket,popen2,time

#
# $Log$
# Revision 1.1  2004/02/20 09:24:34  kpalin
# Version that sort of takes care of the file serving by it self.
# Doesn't really work.
#
#

try:
    from SimpleXMLRPCServer import CGIXMLRPCRequestHandler
except ImportError:
    from SimpleXMLRPCServer233 import CGIXMLRPCRequestHandler

    
hostname=socket.gethostname()
portno=8001


errstrm=open("mabsCom.err","a")
try:
    import socket,popen2,os


    from xmlrpclib import ServerProxy


    servName=socket.gethostname()
    servPort=int(open("mabsCommandServer.portno").read())


    server=ServerProxy("http://%s:%d/"%(servName,servPort))

    def nextFasta(clientID):
        try:
            filen=server.nextFasta(clientID)
            print >> errstrm, "Serving",filen,"to",os.environ["REMOTE_HOST"],"(",clientID,") at", time.asctime()
            errstrm.flush()
            return  filen
        
        except Exception,e:
            print >> errstrm,e
            return "xTryAgain"


    handler = CGIXMLRPCRequestHandler()
    handler.register_function(nextFasta)
    handler.handle_request()
except Exception,e:
    print >> errstrm,e
