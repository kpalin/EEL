"""Takes care of the output formating.

Also gzip related."""

from time import localtime
from cStringIO import StringIO





#
# $Log$
# Revision 1.9  2004/02/26 11:28:22  kpalin
# Changed the output routines to handle the new output from alignment
# bestAlignments
#
# Revision 1.8  2004/02/23 12:23:52  kpalin
# Updates for per gene orthologous runs. Maybe litle multiple alignment.
#
# Revision 1.7  2004/01/23 12:52:00  kpalin
# Some, no doubt vital, modifications
#
# Revision 1.6  2004/01/14 10:05:58  kpalin
# Generated documentation
#
# Revision 1.5  2004/01/09 10:07:00  kpalin
# Output in Anchor format.
# GFF format gives score changes, not increasing sequence of scores.
#
# probably much more
#
# Revision 1.4  2003/12/30 11:20:45  kpalin
# In interface.py Reverse the previous changes and allow passing of
# gziped files to c++ align extension
#
# Revision 1.3  2003/12/29 12:43:18  kpalin
# Interface class repaired to enable alignment from gzip:ed temporary files.
#
# Ilmeisesti jotain uutta. En tiedä mitä.
#
# Revision 1.2  2003/12/12 12:43:54  kpalin
# Lisättiin loki ja gff2aplot ohjelmalle sopiva tuloste.
# Nyt saa kivoja 2D postscript kuvia helposti.
#
#




try:
    from gzip import GzipFile
except ImportError:
    print "No gzip available."

def savematch(data, filename=''):
    """data must have the following format:
    dictionary from Matrix to Sequence to Index to Score"""
    #Maybe one should add a security policy for allowed filenames.
    #e.g. do not allow '/' in filename.
    if filename=='':
        a=localtime()
        filename='mabs_'+str(a.tm_year)+'_'+str(a.tm_mon)+'_'+str(a.tm_mday)+'_'+str(a.tm_hour)+'_'+str(a.tm_min)+'.gff'
    try:
        if filename[-3:]==".gz":
            try:
                F=GzipFile(filename,"w")
            except NameError:
                filename=filename[:-3]
                F=open(filename,'w')
        else:
                F=open(filename,'w')
        for Matr in data.keys():
            for Seq in data[Matr].keys():
                for Pos,Strand in data[Matr][Seq].keys():
                    F.write("%s\tmabs\t%s\t%d\t%d\t%f\t%s\t.\n"%(Seq,Matr.getName(),Pos,Pos+len(Matr)-1,data[Matr][Seq][(Pos,Strand)],Strand))
        F.close()
        return filename
        
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        return ''

def showmatch(data):
    """data must have the following format:
    dictionary from Matrix to Sequence to Index to Score"""
    for Matr in data.keys():
        for Seq in data[Matr].keys():
            for Pos,Strand in data[Matr][Seq].keys():
                print "%s\tmabs\t%s\t%d\t%d\t%f\t%s\t."%(Seq,Matr.getName(),Pos,Pos+len(Matr)-1,data[Matr][Seq][(Pos,Strand)],Strand)

import types

def get(data):
    """Returns the data as a string formated as GFF.

    data must have the following format:
    dictionary from Matrix to Sequence to Index to Score"""
    #output=''
    def flatten(mlist):
        olist=[]
        for i in mlist:
            if type(i)==types.ListType:
                olist.extend(flatten(i))
            else:
                olist.append(i)
        return olist
    output="".join(flatten(map(lambda Mat:map(lambda Seq:map(lambda ((Ind,strand),score):\
        ("%s\tmabs\t%s\t%d\t%d\t%f\t%s\t.\n"%(Seq,Mat.getName(),Ind,Ind+len(Mat)-1,score,strand)),\
        data[Mat][Seq].items()), data[Mat].keys()),data.keys())))
    return output




def savealign(alignment, filename,mode="w"):
    "saves the alignment"
    try:
        F=open(filename,mode)
        F.write(alignment)
        F.close()
        return filename
        
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        return ''







def formatalignCHAOS(alignment):
    "Formats the alignment as anchoring file for DIALIGN2 (a la CHAOS) Not quite"
    if not alignment:
        return "No alignment\n"

    outStrIO=StringIO()
    xname,yname=alignment.x_name,alignment.y_name

    for i,goodAlign in zip(range(1,len(alignment.bestAlignments)+1),alignment.bestAlignments):
        if len(goodAlign)==0:
            continue
        prevScore=0.0
        #        for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
        for as in goodAlign: # as == Aligned site
            prevScore,score=as.score,as.score-prevScore
            outStrIO.write("\t".join(map(str,[xname,yname,as.beginX,as.beginY,as.endX-as.beginX,score]))+"\n")

    return outStrIO.getvalue()


def formatalignGFF(alignment):
    "Formats the alignment for GFF file"
    if not alignment:
        return "No alignment\n"

    outStrIO=StringIO()
    GFFformat="%s\tmalign\t%s\t%d\t%d\t%4.2f\t%s\t.\tCM %d\tMW %g\n" # CM= Cis Module, MW = Matrix weight
    # This format should be OK for gff2aplot 2.0
    GFFalignFormat='%s\tmalign\t%s\t%d\t%d\t%4.2f\t%s\t.\t\tTarget "%s";\tStart %d;\tEnd %d;\tStrand .;\tFrame .;\tCM %d;\n'

    outStrIO.write("### lambda=%f mu=%f nu=%f xi=%f Nucleotides per rotation=%f time=%g\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation,alignment.secs_to_align))
    xname,yname=alignment.x_name,alignment.y_name

    for i,goodAlign in zip(range(1,len(alignment.bestAlignments)+1),alignment.bestAlignments):
        if len(goodAlign)==0:
            continue
        outStrIO.write(GFFalignFormat%(xname,"CisModule",goodAlign[0].beginX,goodAlign[-1].endX,goodAlign[-1].score,".",yname,goodAlign[0].beginY,goodAlign[-1].endY,i))
        outStrIO.write(GFFalignFormat%(yname,"CisModule",goodAlign[0].beginY,goodAlign[-1].endY,goodAlign[-1].score,".",xname,goodAlign[0].beginX,goodAlign[-1].endX,i))

        prevScore=0.0
        for as in goodAlign:
            prevScore,score=as.score,as.score-prevScore
            outStrIO.write(GFFformat%(xname,as.motif,as.beginX,as.endX,score,as.strand,i,as.siteScoreX))
            outStrIO.write(GFFformat%(yname,as.motif,as.beginY,as.endY,score,as.strand,i,as.siteScoreY))

    return outStrIO.getvalue()


def formatalign(alignment,seq=None):
    "Formats the alignment for human use"
    if not alignment:
        return "No alignment\n"

    outStrIO=StringIO()
    
    try:
        from editdist import alignSeq
    except ImportError:
        print "Using subb alignment"
        def alignSeq(xseq,yseq):
            "Stubb for alignment. Does nothing really"
            l=max(len(xseq),len(yseq))
            return (0,xseq.ljust(l).replace(" ","-"),yseq.rjust(l).replace(" ","-"))


            
        
    def formatAlnSeq(xaln,yaln,xname,yname,xstart,ystart,linelen=60):
        """Formats the DNA sequence alignment output"""
        outstr="Sequence 1: %s\nSequence 2: %s\n\n"%(xname,yname)
        n=len(xaln)
        xpos,ypos=xstart+1,ystart+1
        for i in range(0,n,linelen):
            xline,yline=xaln[i:i+linelen],yaln[i:i+linelen]
            outstr+="%7d : %s\n%7d : %s\n\n"%(xpos,xline,ypos,yline)
            xpos+=len(xline)-xline.count("-")
            ypos+=len(yline)-yline.count("-")
        return outstr

    outStrIO.write("### lambda=%g mu=%g nu=%g xi=%g Nucleotides per rotation=%g\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation))
    outStrIO.write("### D[%s][%s]\nNote! First nucleotide at position 1 (one) and binding site at zero!\n"%(alignment.x_name,alignment.y_name))

    xname,yname=alignment.x_name,alignment.y_name
    xe,ye=0,0
    if seq and not (seq.has_key(xname) and seq.has_key(yname)):
        seq=None

    # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
    for i,goodAlign in zip(range(1,len(alignment.bestAlignments)+1),alignment.bestAlignments):
        if len(goodAlign)==0:
            continue
        if seq:
            xstart=max(goodAlign[0].beginX-10,0)
            ystart=max(goodAlign[0].beginY-10,0)
            xend=min(goodAlign[-1].endX+10,len(seq[xname]))
            yend=min(goodAlign[-1].endY+10,len(seq[yname]))

            xseq=seq[xname][xstart:xend]
            yseq=seq[yname][ystart:yend]

            xaln=""
            yaln=""
            xadded=xstart
            yadded=ystart


        outStrIO.write("\n### Alignment No %d ###\n"%(i,))
        #for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
        for as in goodAlign:
            outStrIO.write("D[%d][%d]=%.2f %s (%d,%d) <=> (%d,%d) %s\n"%(as.seqX,as.seqY,as.score,as.motif,as.beginX,as.endX,as.beginY,as.endY,as.strand))
            if seq:
                y2add=yseq[yadded-ystart:as.beginY-1-ystart].lower()
                x2add=xseq[xadded-xstart:as.beginX-1-xstart].lower()
                #alnFmt="%%s%%-%ds%%s"%(max(len(y2add),len(x2add)))
                alnFmt="%s%s%s"
                distYX,y2add,x2add=alignSeq(y2add,x2add)
                yaln=alnFmt%(yaln,y2add,yseq[as.beginY-1-ystart:as.endY-ystart].upper())
                xaln=alnFmt%(xaln,x2add,xseq[as.beginX-1-xstart:as.endX-xstart].upper())
                xadded,yadded=as.endX,as.endY

        if seq:
            distYX,y2add,x2add=alignSeq(yseq[yadded-ystart:yend-ystart].lower(),xseq[xadded-xstart:xend-xstart].lower())
            yaln+=y2add
            xaln+=x2add

            outStrIO.write("\n"+formatAlnSeq(xaln.replace(" ","-"),yaln.replace(" ","-"),xname,yname,xstart,ystart))
            #outstr+="%s\n%s\n"%(xaln.replace(" ","-"),yaln.replace(" ","-"))

    outStrIO.write("### Alignment took %.1f CPU seconds.\n"%(alignment.secs_to_align))


    return outStrIO.getvalue()

        
