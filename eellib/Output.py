from time import localtime
from cStringIO import StringIO





#
# $Log$
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
##                    F.write(Seq                    +'\t'+  #seqname
##                            'mabs'               +'\t'+  #Add program name here
##                            M.getName()          +'\t'+  #feature
##                            str(Pos)               +'\t'+  #start
##                            str(Pos+len(Matr)-1)      +'\t'+  #end
##                            str(data[Matr][Seq][(Pos,Strand)])+'\t'+  #score
##                            str(Strand)+'\t'+  #strand
##                            '.'                  +'\n')  #frame
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
##                print Seq                    +'\t'+  #seqname
##                        'mabs'               +'\t'+  #Add program name here
##                        M.getName()          +'\t'+  #feature
##                        str(Pos)               +'\t'+  #start
##                        str(Pos+len(Matr)-1)      +'\t'+  #end
##                        str(data[Matr][Seq][(Pos,Strand)])+'\t'+  #score
##                        str(Strand)+'\t'+  #strand
##                            '.'                  +'\n') #frame

##    for M in data.keys():
##        for S in data[M].keys():
##            for I in data[M][S].keys():
##                print(S                    +'\t'+  #seqname
##                      'mabs'               +'\t'+  #Add Program Name here
##                      M.getName()          +'\t'+  #feature
##                      str(I)               +'\t'+  #start
##                      str(I+len(M)-1)      +'\t'+  #end
##                      str(data[M][S][I][0])+'\t'+  #score
##                      str(data[M][S][I][1])+'\t'+  #strand
##                      '.'                  )  #frame

import types

def get(data):
    """data must have the following format:
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
##    for M in data.keys():
##        for S in data[M].keys():
##            for I in data[M][S].keys():
##                output+=(S                    +'\t'+  #seqname
##                         'mabs'               +'\t'+  #Add Program Name here
##                         M.getName()          +'\t'+  #feature
##                         str(I)               +'\t'+  #start
##                         str(I+len(M)-1)      +'\t'+  #end
##                         str(data[M][S][I][0])+'\t'+  #score
##                         str(data[M][S][I][1])+'\t'+  #strand
##                         '.'                  +'\n')  #frame
    return output




def savealign(alignment, filename=''):
    "saves the alignment"
    try:
        #print filename
        F=open(filename,'w')
        F.write(alignment)
        F.close()
        return filename
        
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        return ''



def formatalignGFF(alignment):
    "Formats the alignment for GFF file"
    if not alignment:
        return "No alignment\n"

    outStrIO=StringIO()
    GFFformat="%s\tmalign\t%s\t%d\t%d\t%4.2f\t%s\t.\tCM %d\n"
    # This format should be OK for gff2aplot 2.0
    GFFalignFormat='%s\tmalign\t%s\t%d\t%d\t%4.2f\t%s\t.\t\tTarget "%s";\tStart %d;\tEnd %d;\tStrand +;\tFrame .;\tCM %d;\n'

    outStrIO.write("### lambda=%f mu=%f nu=%f xi=%f Nucleotides per rotation=%f\n"%(alignment.Lambda,alignment.Mu,alignment.Nu,alignment.Xi,alignment.nuc_per_rotation))
    xname,yname=alignment.x_name,alignment.y_name

    for i,goodAlign in zip(range(1,len(alignment.bestAlignments)+1),alignment.bestAlignments):
        if len(goodAlign)==0:
            continue
        outStrIO.write(GFFalignFormat%(xname,"CisModule",goodAlign[0][4][0],goodAlign[-1][4][1],goodAlign[-1][2],".",yname,goodAlign[0][5][0],goodAlign[-1][5][1],i))
        outStrIO.write(GFFalignFormat%(yname,"CisModule",goodAlign[0][5][0],goodAlign[-1][5][1],goodAlign[-1][2],".",xname,goodAlign[0][4][0],goodAlign[-1][4][1],i))
        for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
            outStrIO.write(GFFformat%(xname,motif,xcoord[0],xcoord[1],score,strand,i))
            outStrIO.write(GFFformat%(yname,motif,ycoord[0],ycoord[1],score,strand,i))

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
            xstart=max(goodAlign[0][4][0]-10,0)
            ystart=max(goodAlign[0][5][0]-10,0)
            xend=min(goodAlign[-1][4][1]+10,len(seq[xname]))
            yend=min(goodAlign[-1][5][1]+10,len(seq[yname]))

            xseq=seq[xname][xstart:xend]
            yseq=seq[yname][ystart:yend]

            xaln=""
            yaln=""
            xadded=xstart
            yadded=ystart


        outStrIO.write("\n### Alignment No %d ###\n"%(i,))
        for (x,y,score,motif,xcoord,ycoord,strand) in goodAlign:
            outStrIO.write("D[%d][%d]=%.2f %s (%d,%d) <=> (%d,%d) %s\n"%(x,y,score,motif,xcoord[0],xcoord[1],ycoord[0],ycoord[1],strand))
            if seq:
                y2add=yseq[yadded-ystart:ycoord[0]-1-ystart].lower()
                x2add=xseq[xadded-xstart:xcoord[0]-1-xstart].lower()
                #alnFmt="%%s%%-%ds%%s"%(max(len(y2add),len(x2add)))
                alnFmt="%s%s%s"
                distYX,y2add,x2add=alignSeq(y2add,x2add)
                yaln=alnFmt%(yaln,y2add,yseq[ycoord[0]-1-ystart:ycoord[1]-ystart].upper())
                xaln=alnFmt%(xaln,x2add,xseq[xcoord[0]-1-xstart:xcoord[1]-xstart].upper())
                xadded,yadded=xcoord[1],ycoord[1]

        if seq:
            distYX,y2add,x2add=alignSeq(yseq[yadded-ystart:yend-ystart].lower(),xseq[xadded-xstart:xend-xstart].lower())
            yaln+=y2add
            xaln+=x2add

            outStrIO.write("\n"+formatAlnSeq(xaln.replace(" ","-"),yaln.replace(" ","-"),xname,yname,xstart,ystart))
            #outstr+="%s\n%s\n"%(xaln.replace(" ","-"),yaln.replace(" ","-"))

    outStrIO.write("### Alignment took %.1f CPU seconds.\n"%(alignment.secs_to_align))


    return outStrIO.getvalue()

        
