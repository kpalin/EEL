#!/usr/bin/python2.2
import socket,popen2,time

#
# $Log$
# Revision 1.2  2004/02/20 10:57:26  kpalin
# For use in apache web server so that the distribution
# can go through firewalls.
#
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
    (servPort,servPid,servComm)=eval(open("mabsCommandServer.portno").read())


    server=ServerProxy("http://%s:%d/"%(servName,servPort))

    def nextFasta(clientID):
        try:
            filen=server.nextFasta(clientID)
            print >> errstrm, "Serving",filen,"to",os.environ["REMOTE_HOST"],"(",clientID,") at", time.asctime()
            errstrm.flush()
            return  filen
        
        except Exception,e:
            print >> errstrm,e
            if e[0]==111 and os.system("/bin/ps -no-headers %d >/dev/null"%(servPid)) != 0:
                print >> errstrm, "Looks like server has died. Trying to start it."
                os.system("python2.2 %s 2>/dev/null >/dev/null </dev/null &"%(servComm))
                
            return "xTryAgain"


    handler = CGIXMLRPCRequestHandler()
    handler.register_function(nextFasta)
    handler.handle_request()
except Exception,e:
    print >> errstrm,e
