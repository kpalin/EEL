#!/usr/bin/python2.2
from CGIXMLRPCServer import SimpleXMLRPCServer233
import socket,popen2

#
# $Log$
#



hostname=socket.gethostname()
portno=8001

FileList=eval(open("serverFiles.py").read())

FileList.reverse()

def fileExists(fname):
    try:
        os.stat(fname)
        return 1
    except OSError:
        return 0

import sys,os
def nextFasta():
    global FileList
    try:
        while 1:
            filen=FileList.pop()
            if (not fileExists("/home/kpalin/fs/tyot/mabs/synteny/output/good/%s.align.gff.gz"%filen)) and \
               (not fileExists("/home/kpalin/fs/tyot/mabs/synteny/output/%s.align.gff.gz"%filen)):
                break
            else:
                print >> sys.stderr, "Skipping",filen
        print >> sys.stderr, filen
        return filen
    except IndexError:
        return "xAllDone"
    

handler = CGIXMLRPCRequestHandler()
handler.register_function(nextFasta)
handler.handle_request()
