#!/usr/bin/python2.2
from SimpleXMLRPCServer import SimpleXMLRPCServer
import socket,popen2
import sys,os

#
#  $Log$
#  Revision 1.2  2004/02/23 12:22:57  kpalin
#  Updates for per gene orthologous runs.
#
#  Revision 1.1  2004/02/20 10:56:05  kpalin
#  Serves the list of files given as parameter.
#  Keeps some log who is running what.
#
#



hostname=socket.gethostname()
portno=8001


try:
    FileListName=sys.argv[1]
except IndexError:
    print "Usage: python mabsCommandServer.py listOfFilesToServe.py"
runListName="%sOnRun.py"%(FileListName)

FileList=eval(open(FileListName).read())

def fileExists(dirs,fname):
    for d in dirs:
        pname="%s/%s"%(d,fname)
        try:
            os.stat(pname)
            print >> sys.stderr, "match",pname
            return 1
        except OSError:
            pass
    return 0

try:
    runningFiles=eval(open(runListName).read())
except IOError:
    runningFiles={}

supposedlyDone={}

def nextFasta(clientID):
    global FileList,runListName,supposedlyDone
    dirs=["/fs/home/kpalin/tyot/mabs/synteny/output/good/",
          "/fs/home/kpalin/tyot/mabs/synteny/output/",
          "/fs/home/kpalin/tyot/mabs/synteny/output/good/ws/",
          "/fs/home/kpalin/tyot/mabs/synteny/perGene/output/",
          "/fs/home/kpalin/tyot/mabs/synteny/perGene/output/ws/output/"]
    try:
        
        for i in range(len(FileList)):
            filen=FileList[i]
            if  fileExists(dirs,filen+".align.gff.gz"):
                FileList[i]=None
            elif (not filen in runningFiles.values()):
                break
            print >> sys.stderr, "Skipping",filen
        print >> sys.stderr, "Serving",filen

        FileList=filter(lambda x:x,FileList)
        if runningFiles.has_key(clientID):
            i=FileList.index(runningFiles[clientID])
            print >> sys.stderr, "Looks like done",FileList[i]
            del FileList[i]
        runningFiles[clientID]=filen

        ## To be sure it's OK.
        open(FileListName,"w").write(repr(FileList))
        open(runListName,"w").write(repr(runningFiles))
        return filen
    except Exception,e:
        import traceback
        traceback.print_exc(file=sys.stderr)
        return "xAllDone"
    


OK=0
while not OK:
    try:
        print "Starting server at http://%s:%d"%(hostname,portno)
        server=SimpleXMLRPCServer((hostname,portno))
        OK=1
    except socket.error,e:
        portno=portno+1
        print e[1]


open("mabsCommandServer.portno","w").write(repr((portno,os.getpid()," ".join(sys.argv))))

server.register_function(nextFasta)
server.serve_forever()
