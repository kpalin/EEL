# -*- coding: UTF-8 -*-

"""Takes care of the output formating.

Also gzip related."""

from time import localtime
from cStringIO import StringIO

from eellib import alignedCols
import math


#
# $Log$
# Revision 1.32  2010-05-14 11:10:42  sksipila
# Added shortMultipleAlign, multiFromPairwise and treeMultipleAlign
#
# Revision 1.31  2008/03/12 12:03:03  kpalin
# Fixed "documentation"
#
# Revision 1.30  2008/02/29 10:08:23  kpalin
# Cut motif name from the end.
#
# Revision 1.29  2008/01/02 08:27:14  kpalin
# Added sequence names and descriptions to GFF output file.
#
# Revision 1.28  2006/11/13 13:03:01  kpalin
# Added module E-score RMSE
#
# Revision 1.27  2006/11/13 12:33:41  kpalin
# Added 50mer tags to the GFF output giving the location of the
# module. The tags start from the beginning and the end of the module.
#
# Added Escore reporting.
#
# Revision 1.26  2006/08/31 10:08:36  kpalin
# More Debug output and TF distance calculations.
#
# Revision 1.25  2006/08/14 09:46:06  kpalin
# Added more comments.
#
# Revision 1.24  2006/02/08 07:44:46  kpalin
# Fixed faulty asserts.
#
# Revision 1.23  2005/11/25 12:09:55  kpalin
# Human readable format for multiple alignment.
#
# Revision 1.22  2005/11/24 13:29:13  kpalin
# Added some debugging code for GFF alignment output.
#
# Revision 1.21  2005/05/19 07:49:35  kpalin
# Merged Waterman-Eggert style suboptimal alignments and
# SNP matching.
#
# Revision 1.20.2.2  2005/05/10 13:12:14  kpalin
# Added a 3rd row for the alignment noting the aligned motif.
#
# Revision 1.20.2.1  2005/04/12 09:12:10  kpalin
# Formatalign outputs also the description of the sequences, when
# available.
#
# Revision 1.20  2005/03/03 09:04:59  kpalin
# Even better output for SNP matching.
#
# Revision 1.19  2005/02/24 11:37:43  kpalin
# Site annotations.
#
# Revision 1.18  2005/02/22 11:08:56  kpalin
# Initial working version for outputting results from
# SNP scanning matrix matcher.
#
# Revision 1.17  2005/01/13 13:16:42  kpalin
# Moved the requesting of sequences to be aligned to Python side
# of stuff. Much better.
#
# Revision 1.16  2005/01/12 13:34:55  kpalin
# Added Tkinter/Tix Graphical user interface and command -no-gui to
# avoid it.
#
# Revision 1.15  2005/01/07 13:41:25  kpalin
# Works with py2exe. (windows executables)
#
# Revision 1.14  2004/12/17 12:21:26  kpalin
# Changed the matrix dictionary key ordering and made the GFF output
# routines routine.
#
# Revision 1.13  2004/12/14 13:08:05  kpalin
#
# Name change from MABS to EEL (Enhancer Element Locator / Monty Python pun
# "My hovercraft is full of EELs" )
#
# Revision 1.12  2004/10/21 12:45:43  kpalin
# Added possibility to gzip also the alignment files.
#
# Revision 1.11  2004/07/30 12:08:54  kpalin
# Changed the GFF formatting to use alnColumn objects. (Needed
# for multiple alignment)
#
# Revision 1.10  2004/03/03 09:26:59  kpalin
# Few possible bugs corrected.
#
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



try: # Workaround for python 2.2
    enumerate([2,3,4])
except NameError:
    enumerate=lambda x:zip(range(len(x)),x)


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
        filename='eel_'+str(a.tm_year)+'_'+str(a.tm_mon)+'_'+str(a.tm_mday)+'_'+str(a.tm_hour)+'_'+str(a.tm_min)+'.gff'
    try:
        if filename[-3:]==".gz":
            try:
                F=GzipFile(filename,"w")
            except NameError:
                filename=filename[:-3]
                F=open(filename,'w')
        else:
                F=open(filename,'w')

## This is in wrong format Seq and Matr are reversed.
##        for Matr in data.keys():
##            for Seq in data[Matr].keys():
##                for Pos,Strand in data[Matr][Seq].keys():
##                    F.write("%s\teel\t%s\t%d\t%d\t%f\t%s\t.\n"%(Seq,Matr.getName(),Pos,Pos+len(Matr)-1,data[Matr][Seq][(Pos,Strand)],Strand))
        F.write(get(data))
        F.close()
        return filename
        
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        return ''

def showmatch(data):
    """data must have the following format:
    dictionary from Matrix to Sequence to Index to Score"""
##    for Matr in data.keys():
##        for Seq in data[Matr].keys():
##            for Pos,Strand in data[Matr][Seq].keys():
##                print "%s\teel\t%s\t%d\t%d\t%f\t%s\t."%(Seq,Matr.getName(),Pos,Pos+len(Matr)-1,data[Matr][Seq][(Pos,Strand)],Strand)

    if type(data)==type(""):
        print data
    else:
        print get(data)

import types

def get(data):
    """Returns the data as a string formated as GFF.

    data must have the following format:
    dictionary from Sequence to Matrix to (position,strand,((ambig,allele,snpPos,scoreDif))) to (Score,altScore)"""
    #output=''
    # data = {'Sequence':{Matrix:{(position,strand,((ambig,allele,snpPos,scoreDif))):(Score,altScore)}}}

    def flatten(mlist):
        olist=[]
        for i in mlist:
            if type(i)==types.ListType:
                olist.extend(flatten(i))
            else:
                olist.append(i)
        return olist
    output="".join(flatten(map( \
        lambda Seq:map(lambda Mat:map(lambda ((Ind,strand,snps),(score,altScore)):\
        ("%s\teel\t%s\t%d\t%d\t%f\t%s\t.\t%s\n"%(\
        Seq,Mat.getName(),Ind,Ind+len(Mat)-1,score,strand,\
        "\t".join(["%s%s%d %g"%(snp[0],snp[1],snp[2],snp[3]) for snp in snps]))),\
        data[Seq][Mat].items()), data[Seq].keys()),data.keys())))
    return output




def savealign(alignment, filename,mode="w"):
    "saves the alignment"
    try:
        if filename[-3:]==".gz":
            try:
                F=GzipFile(filename,mode)
            except NameError:
                filename=filename[:-3]
                F=open(filename,mode)
        else:
            F=open(filename,mode)
        F.write(alignment)
        F.close()
        return filename
        
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        return ''





def moduleData(AlnList):
    """Returns list of information about the alignment module.
    Input is a list of alignedCols.alnColumn objects.
    Output is a list of [sequencecode,begin,end,score] lists"""
    byCodeMap={}
    for col in AlnList:
        for seq,(begin,end) in zip(col.seqID,col.beginEnd):
            if byCodeMap.has_key(seq):
                byCodeMap[seq][1]=min(begin,byCodeMap[seq][1])
                byCodeMap[seq][2]=max(end,byCodeMap[seq][2])
                byCodeMap[seq][3]=max(col.score,byCodeMap[seq][3])
            else:
                byCodeMap[seq]=[seq,begin,end,col.score]
    outp=byCodeMap.values()
    outp.sort()
    return outp


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


def formatMultiAlignGFF(alignment,seqData=None,tagLength=50):
    "Formats the multiple alignemnt for GFF file"
    if not alignment:
        return "No alignment\n"
    
    outStrIO=StringIO()
    
    # Header
    NamesAndLengths=""
    #NamesAndLengths=["%s=%d"%(x,y) for (x,y) in zip(alignment.names,alignment.lengths) ]
    outStrIO.write("### lambda=%f mu=%f nu=%f xi=%f Nucleotides per rotation=%f time=%g %s\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation,alignment.secs_to_align,NamesAndLengths))
    
    for sName in alignment.names:
        try:
            sDesc=seqData.describe(sName).split("\n")
            sDesc="\n#".join(sDesc)
        except (KeyError,AttributeError):
            sDesc=""
        outStrIO.write("#%s\t%s\n"%(sName,sDesc))
    
    # CM= Cis Module, MW = Matrix weight, COL = Column code on module
    GFFformat="%s\tmalign\t%s\t%d\t%d\t%4.2f\t%s\t.\tCM %d\tMW %g\tCOL %d\t%s\n" 
    
    # the last %s is for inserting parameters for gff2aplot etc.
    GFFalignFormat='%s\tmalign\t%s\t%d\t%d\t%4.2f\t.\t.\t%sCM %d\n'
    # For each module
    for i,goodAlign in zip(range(1,len(alignment.bestAlignments)+1),alignment.bestAlignments):
	if len(goodAlign)==0:
            continue
        DEBUGprevpos={}
        
        # For each sequence participating in the module
        modData=moduleData(goodAlign)
        for j in range(len(modData)):
            seq,begin,end,score=modData[j]
            InsertInfo=""
            if len(modData)==2:
                targetID=1-abs(j)
                InsertInfo='\tTarget "%s";\tStart %d;\tEnd %d;\tStrand .;\tFrame .;\t%s'%(alignment.names[modData[targetID][0]],\
                                                                                       modData[targetID][1],modData[targetID][2],InsertInfo)
            
            
            if alignment.Rsquared>0.1:
                InsertInfo="%sEscore %g\t"%(InsertInfo,math.exp(alignment.alpha+alignment.beta*score))
            
            if seqData:
                moduleSeq=seqData[alignment.names[seq]][begin:end]
                InsertInfo="%sTag5 %s\tTag3 %s\t"%(InsertInfo,moduleSeq[:tagLength],moduleSeq[-tagLength:])
            
            else:
                InsertInfo=""
            outStrIO.write(GFFalignFormat%(alignment.names[seq], \
                                           "CisModule",\
                                           begin,end,\
                                           score,InsertInfo,i))
        
        # For each column on the module.
        prevScore=0.0
        for as,colCode in zip(goodAlign,range(len(goodAlign))):
            prevScore,score=as.score,as.score-prevScore
            for seq,(begin,end),siteScore,annot in zip(as.seqID,as.beginEnd,as.siteScore,as.annotation):
                if type(siteScore)==type(1.1):
                    outStrIO.write(GFFformat%(alignment.names[seq], \
                                              as.motif, \
                                              begin,end,score,as.strand,\
                                              i,siteScore,colCode,annot.strip()))
                    try:
                        assert(DEBUGprevpos.get(seq,(0,0))[1]<begin)
                        assert(1<=begin-DEBUGprevpos.get(seq,(0,begin-1))[1]<=1000)
                    except AssertionError:
                        print DEBUGprevpos.get(seq,(0,0))[1],"<",begin
                        print "0<=",begin,"-",DEBUGprevpos.get(seq,(0,begin-1))[1],"<=",1000
                        print alignment.names[seq],"COL %d"%(colCode),DEBUGprevpos.get(seq,(0,0))[1],"<",end
#                        print outStrIO.getvalue()
                        raise
                    DEBUGprevpos[seq]=(begin,end)
    
    return outStrIO.getvalue()


def formatalignGFF(alignment):
    "Formats the alignment for GFF file"
    return formatMultiAlignGFF(alignment)



def formatalign2D(alignment,seq=None):
    "Formats the alignment for human use"
    if not alignment:
        return "No alignment\n"

    outStrIO=StringIO()
    
    try:
        from editdist import alignSeq
    except ImportError:
        print "Using subb alignment"
        def alignSeq(xseq,yseq,*rest):
            "Stubb for alignment. Does nothing really"
            l=max(len(xseq),len(yseq))
            return (0,xseq.ljust(l).replace(" ","-"),yseq.rjust(l).replace(" ","-"))


            
        
    def formatAlnSeq(xaln,yaln,xname,yname,xstart,ystart,motifAln="",linelen=60):
        """Formats the DNA sequence alignment output"""
        outstr="Sequence 1: %s\nSequence 2: %s\n\n"%(xname,yname)
        n=len(xaln)
        xpos,ypos=xstart+1,ystart+1
        for i in range(0,n,linelen):
            xline,yline,mline=xaln[i:i+linelen],yaln[i:i+linelen],motifAln[i:i+linelen]
            outstr+="%7d : %s\n%7d : %s\n          %s\n\n"%(xpos,xline,ypos,yline,mline)
            xpos+=len(xline)-xline.count("-")
            ypos+=len(yline)-yline.count("-")
##        n=len(xaln)
##        xpos,ypos=xstart+1,ystart+1
##        for i in range(0,n,linelen):
##            xline,yline=xaln[i:i+linelen],yaln[i:i+linelen]
##            outstr+="%7d : %s\n%7d : %s\n\n"%(xpos,xline,ypos,yline)
##            xpos+=len(xline)-xline.count("-")
##            ypos+=len(yline)-yline.count("-")
        return outstr

    outStrIO.write("### lambda=%g mu=%g nu=%g xi=%g Nucleotides per rotation=%g\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation))
    outStrIO.write("### D[%s][%s]\nNote! First nucleotide at position 1 (one) and binding site at zero!\n"%(alignment.x_name,alignment.y_name))

    xname,yname=alignment.x_name,alignment.y_name
    xe,ye=0,0
    if seq and not (seq.has_key(xname) and seq.has_key(yname)):
        seq=None

    if seq:
        outStrIO.write("Sequence %s:\n%s\n\nSequence %s:\n%s\n"%(xname,seq.describe(xname),yname,seq.describe(yname)))
        
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
            maln=""
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
                #alnFmt="%s%s%s"
                siteLen=as.endX-as.beginX+1
                alnFmt="%%s%%s%%-%ds"%(siteLen)

                assert(len(xseq[as.beginX-1-xstart:as.endX-xstart])==siteLen)
                distYX,y2add,x2add=alignSeq(y2add,x2add,1,1)
                yaln=alnFmt%(yaln,y2add,yseq[as.beginY-1-ystart:as.endY-ystart].upper())
                xaln=alnFmt%(xaln,x2add,xseq[as.beginX-1-xstart:as.endX-xstart].upper())

                maln=alnFmt%(maln," "*len(y2add),as.motif[-siteLen:])
                xadded,yadded=as.endX,as.endY

        if seq:
            distYX,y2add,x2add=alignSeq(yseq[yadded-ystart:yend-ystart].lower(),xseq[xadded-xstart:xend-xstart].lower())
            yaln+=y2add
            xaln+=x2add

            outStrIO.write("\n"+formatAlnSeq(xaln.replace(" ","-"),yaln.replace(" ","-"),xname,yname,xstart,ystart,maln))
            #outstr+="%s\n%s\n"%(xaln.replace(" ","-"),yaln.replace(" ","-"))

    outStrIO.write("### Alignment took %.1f CPU seconds.\n"%(alignment.secs_to_align))


    return outStrIO.getvalue()

def formatalign(alignment,seq=None):
    "Formats the alignment for human use"
    if not alignment:
        return "No alignment\n"
    
    outStrIO=StringIO()
    
    try:
        assert(len(alignment.names)==2)  # Provoking AssertionError
        from editdist import alignSeq
        def alignSeqs(seqs,*rest):
            distXY,X,Y=alignSeq(seqs[0],seqs[1],*rest)
            return (distXY,[X,Y])
    
    except(ImportError,AssertionError):
        print "Using subb alignment"
    
        def alignSeqs(seqs,*rest):
            "Stubb for alignment. Does nothing really"
            l=max(map(len,seqs))
            return (0,[xseq.ljust(l).replace(" ","-") for xseq in seqs])
    
    def formatAlnSeq(alns,names,starts,motifAln="",linelen=60):
        """Formats the DNA sequence alignment output"""
        outstr="\n".join(["Sequence %d: %s"%(seqI+1,iname) for (seqI,iname) in enumerate(names)])+"\n\n"
        n=len(alns[0])
        poses=[xstart+1 for xstart in starts]
        for i in range(0,n,linelen):
            lines=[xaln[i:i+linelen] for xaln in alns]
            mline=motifAln[i:i+linelen]
            outstr+="\n".join(["%7d : %s"%(xpos,xline) for (xpos,xline) in zip(poses,lines)])
            outstr+="\n          %s\n\n"%(mline)
            poses=[xpos+len(xline)-xline.count("-") for (xpos,xline) in zip(poses,lines)]
        return outstr
    
    outStrIO.write("### lambda=%g mu=%g nu=%g xi=%g Nucleotides per rotation=%g\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation))
    outStrIO.write("### D%s\nNote! First nucleotide at position 1 (one) and binding site at zero!\n"%("".join(["[%s]"%(x) for x in alignment.names])))
    
    #xname,yname=alignment.x_name,alignment.y_name
    xe,ye=0,0
    if seq:  #Must have all sequences available (Could relax this in future)
       for xseq in alignment.names:
            if not seq.has_key(xseq):
                seq=None
                break
    
    if seq:
        outStrIO.write("\n".join(["Sequence %s:\n%s\n"%(xseq,seq.describe(xseq)) for xseq in alignment.names])+"\n\n")
    
    # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
    for alnNo,goodAlign in enumerate(alignment.bestAlignments):
        if len(goodAlign)==0:
            continue
        if seq:
            starts=[max(spos-10,0) for (spos,epos) in goodAlign[0].beginEnd]
            ends=[min(epos+10,len(seq[alignment.names[i]])) for (i,(spos,epos)) in enumerate(goodAlign[-1].beginEnd)]
            
            seqs=[seq[alignment.names[i]][istart:iend] for (i,(istart,iend)) in enumerate(zip(starts,ends))]
            
            alns=[""]*len(seqs)
            maln=""
            addeds=starts[:]
        
        outStrIO.write("\n### Alignment No %d ###\n"%(alnNo+1,))
        
        #for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
        for as in goodAlign:
            coordStr=""
	    seqCoordStr=""
            for i,x,(spos,epos) in zip(range(1,len(as.siteSeqPos)+1),as.siteSeqPos,as.beginEnd):
                if type(x)==type(1):
                    if(i==1):
                        coordStr="".join([coordStr,"[%d]"%(x)])
                        seqCoordStr="".join([seqCoordStr,"(%d,%d)"%(spos,epos)])
                    else:
                        coordStr="".join([coordStr,"[%d]"%(x)])
                        seqCoordStr=" <=> ".join([seqCoordStr,"(%d,%d)"%(spos,epos)])
                else:
                    if(i==1):
                        coordStr="".join([coordStr,"[ ]"])
                        seqCoordStr="".join([seqCoordStr,"( , )"])
                    else:
                        coordStr="".join([coordStr,"[ ]"])
                        seqCoordStr=" <=> ".join([seqCoordStr,"( , )"])

            outStrIO.write("D%s=%.2f %s %s %s\n"%(coordStr,as.score,as.motif,seqCoordStr,as.strand))
            if seq:
                ToAdd=[xseq[xadded-xstart:beginX-1-xstart].lower() for (xseq,xadded,xstart,(beginX,endX)) in zip(seqs,addeds,starts,as.beginEnd)]
                siteLen=as.beginEnd[0][1]-as.beginEnd[0][0]+1
                alnFmt="%%s%%s%%-%ds"%(siteLen)
                
                #assert(len(xseq[as.beginX-1-xstart:as.endX-xstart])==siteLen)
                distYX,ToAdd=alignSeqs(ToAdd,1,1)
                alns=[alnFmt%(xaln,x2add,xseq[beginX-1-xstart:endX-xstart].upper()) for (xaln,x2add,xseq,xstart,(beginX,endX)) in zip(alns,ToAdd,seqs,starts,as.beginEnd)]
                
                maln=alnFmt%(maln," "*len(ToAdd[0]),as.motif[:siteLen])
                addeds=[endX for (beginX,endX) in as.beginEnd]
        
        if seq:
            ToAdd=[seq[xseq][xadded-xstart:xend-xstart].lower() for (xseq,xadded,xstart,xend) in zip(alignment.names,addeds,starts,ends)]
            distYX,ToAdd=alignSeqs(ToAdd)
            alns=["".join(x) for x in zip(alns,ToAdd)]
            
            outStrIO.write("\n"+formatAlnSeq([xaln.replace(" ","-") for xaln in alns],alignment.names,starts,maln))
            #outstr+="%s\n%s\n"%(xaln.replace(" ","-"),yaln.replace(" ","-"))
    
    outStrIO.write("### Alignment took %.1f CPU seconds.\n"%(alignment.secs_to_align))
    return outStrIO.getvalue()

def formatpwbase(alignment):
    "Formats the seqments of pairwise alignment that multiple alignment was based on for human use"
    if not alignment:
        return "No alignment\n"
    
    outStrIO=StringIO()
    
    outStrIO.write("### lambda=%g mu=%g nu=%g xi=%g Nucleotides per rotation=%g\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation))
    outStrIO.write("### D%s\nNote! First nucleotide at position 1 (one) and binding site at zero!\n"%("".join(["[%s]"%(x) for x in alignment.names])))
    
    for groupNo,groupOfAligns in enumerate(alignment.pwBase):
        # goodAlign= [ (x,y,Score,Motif,(startX,endX),(startY,endY),Strand) ]
        for alnNo,align in enumerate(groupOfAligns):
            if len(align)==0:
                continue
            
            outStrIO.write("\n### Alignment No %d of Group %d ###\n"%(alnNo+1,groupNo+1))
            outStrIO.write("[")
            outStrIO.write(alignment.names[align[0].seqID[0]])
            outStrIO.write("][")
            outStrIO.write(alignment.names[align[0].seqID[1]])
            outStrIO.write("]\n")
            
            #for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
            for as in align:
                coordStr=""
                seqCoordStr=""
                match=0
                for i,x,(spos,epos) in zip(range(1,len(as.siteSeqPos)+1),as.siteSeqPos,as.beginEnd):
                    if(i==1):
                        coordStr="".join([coordStr,"[%d]"%(x)])
                        seqCoordStr="".join([seqCoordStr,"(%d,%d)"%(spos,epos)])
                    else:
                        coordStr="".join([coordStr,"[%d]"%(x)])
                        seqCoordStr=" <=> ".join([seqCoordStr,"(%d,%d)"%(spos,epos)])
                
                for gaNo,goodAlign in enumerate(alignment.bestAlignments):
                    if(gaNo/alignment.numofalign==groupNo and match==0):
                        for gas in goodAlign:
                            for (gspos,gepos) in gas.beginEnd:
                                for(spos,epos) in as.beginEnd:
                                    if(spos==gspos and epos==gepos):
                                        match=match+1
                            if(match>1):
                                break
                            else:
                                match=0
                if(match>1):
                    outStrIO.write("!")
                outStrIO.write("D%s=%.2f %s %s %s\n"%(coordStr,as.score,as.motif,seqCoordStr,as.strand))
        outStrIO.write("\n\n")
    
    return outStrIO.getvalue()
