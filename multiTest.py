

# $Log$

files=[#"ENSDARG00000006837.ENSG00000134323.align.gff",
       #"ENSDARG00000006837.ENSMUSG00000037169.align.gff",
       #"ENSDARG00000006837.ENSRNOG00000006308.align.gff",
       #"ENSDARG00000006837.SINFRUG00000128278.align.gff",
       "ENSG00000134323.ENSMUSG00000037169.align.gff",
       "ENSG00000134323.ENSRNOG00000006308.align.gff",
       #"ENSG00000134323.SINFRUG00000128278.align.gff",
       "ENSMUSG00000037169.ENSRNOG00000006308.align.gff",
       #"ENSMUSG00000037169.SINFRUG00000128278.align.gff",
       #"ENSRNOG00000006308.SINFRUG00000128278.align.gff"
       ]



files=["hum_hum_1031_937.abs9.chrIbg.align.gff",
       "hum_mus_1031.abs9.chrIbg.align.gff",
       "mus_hum_1031_937.abs9.chrIbg.align.gff"]

species={"HUMAN|4.119820698-119936581":0,
         "HUMAN|8.97218432-97331876":1,
         "MOUSE|13.63720377-63822119":2}
specCount=0





class AlignedSite:
    def __init__(self,sType,sScore,pos):
        self.type=sType
        self.score=sScore
        self.pos=pos
        self.pairs=[]

    def __getitem__(self,key):
        return self.pos[key]


    def copy(self):
        c=AlignedSite(self.type,self.score,self.pos)
        c.pairs=self.pairs[:]
        return c


    def addPair(self,spec1,spec2):
        pairs=[spec1,spec2]
        pairs.sort()
        self.pairs.append(tuple(pairs))

    def __repr__(self):
        s="AlignedSite(%s,%s,%s)"%(repr(self.type),repr(self.score),repr(self.pos))
        return s

    def __str__(self):
        return self.__repr__()

    def unite(self,other):
        assert(self.type==other.type)
        self.pos.update(other.pos)
        self.pairs.extend(other.pairs)



class CisMod:
    def __init__(self,species=[],code=[]):
        self.species=species
        #self.species.sort()
        if type(code) not in (type([]),type(())):
            code=[code]
        self.cisCode=code
        self.sites=[]
        self.modScore=0.0

    def __str__(self):
        s=StringIO()
        s.write("%s (%s)"%(" <=> ".join(self.species),"<=>".join(map(lambda x:"CM %d"%(x),self.cisCode))))
        for site in self.sites:
            s.write("\n%s %d %s"%(" ".join(map(lambda x:str(site[x]),self.species)),len(site.pairs),site.type))

        return s.getvalue()

    def __len__(self):
        return len(self.sites)


    def addSite(self,pos,sType,sScore):
        self.sites.append(AlignedSite(sType,sScore,pos))
        if sScore>self.modScore: self.modScore=sScore

    def compatible(self,other):
        """Return the name of the common species if the alignments are between
        three different species"""
        compValue=0
        sSpec,cSpec,oSpec="","",""
        for tSpec in self.species:
            try:
                cSpec=other.species.index(tSpec)
                oSpec=other.species[1-cSpec]
                cSpec=other.species[cSpec]
                compValue=compValue+1
            except ValueError:
                sSpec=tSpec
        if compValue==1:
            return (sSpec,cSpec,oSpec)
        else:
            return None

    def intersect(self,other):
        specs=self.compatible(other)
        if not specs:
            return None

        sSpec,cSpec,oSpec=specs



        i,j=0,0
        m=0

        newCis=CisMultiMod([sSpec,cSpec,oSpec],self.cisCode+other.cisCode)
        #s="%s <=> %s <=> %s (CM%d<=>CM%d)\n"%(sSpec,cSpec,oSpec,self.cisCode,other.cisCode)
        s=str(newCis)
        while i<len(self.sites) and j<len(other.sites):
            if self.sites[i][cSpec]<other.sites[j][cSpec]:
                i+=1
            elif self.sites[i][cSpec]>other.sites[j][cSpec]:
                j+=1
            elif self.sites[i].type==other.sites[j].type:
                newSite=self.sites[i].copy()
                newSite.unite(other.sites[j])
                newCis.sites.append(newSite)

                s=s+"%d %d %d %s %s\n"%(self.sites[i][sSpec],self.sites[i][cSpec],other.sites[j][oSpec],str(self.sites[i].type),str(newSite.pairs))
                
                m=m+1
                i+=1
                j+=1                
            else:
                i+=1
                j+=1
        if m>1:
            #print s
            print newCis
            #print "Matches/mismatches: %d/%d"%(2*m,len(self.sites)+len(other.sites)-2*m)
            return newCis
        else:
            return None

from cStringIO import StringIO


class CisMultiMod(CisMod):
    def __init__(self,species=[],CMcodes=[]):
        CisMod.__init__(self,species,CMcodes)


    def addPairAlign(self,other):
        assert(type(other)==type(CisMod()))
        
        spec1,spec2=other.species
        newCis=CisMultiMod(self.species[:],self.cisCode+other.cisCode)

        spair=other.species[:]
        spair.sort()
        try:
            i=self.sites[0].pairs.index(tuple(spair))
            return None
        except ValueError:
            #print self.sites[0].pairs,spair
            pass

        i,j=0,0
        m=0

        while i<len(self.sites) and j<len(other.sites):
            if self.sites[i][spec1]<other.sites[j][spec1] or self.sites[i][spec2]<other.sites[j][spec2]:
                i+=1
            elif self.sites[i][spec1]>other.sites[j][spec1] or self.sites[i][spec2]>other.sites[j][spec2]:
                j+=1
            elif self.sites[i].type==other.sites[j].type:
                newSite=self.sites[i].copy()
                newSite.unite(other.sites[j])
                newCis.sites.append(newSite)

                m=m+1
                i+=1
                j+=1                
            else:
                i+=1
                j+=1
        if m>1:
            #print newCis
            #print "Matches/mismatches: %d/%d"%(2*m,len(self.sites)+len(other.sites)-2*m)
            return newCis
        else:
            return None













import re

CMpat=re.compile(r"CM (\d+)")
aligns={}
CisMods={}
nodesInCis={}

cisModules=[]
for file in files:
    handle=open(file)
    print file
    parts=file.split(".")
    if not species.has_key(parts[0]):
        species[parts[0]]=specCount
        specCount+=1
    if not species.has_key(parts[1]):
        species[parts[1]]=specCount
        specCount+=1
    handle.readline()
    line=handle.readline()  # First CisModule
    parts=line.split("\t")
    while len(line)>0:
        spec1=species[parts[0]]
        spec1str=parts[0]
        
        CMcode=CMpat.search(parts[-1])
        CMcode=int(CMcode.group(1))
        #print CMcode
        line=handle.readline() #Second CisModule
        parts=line.split("\t") 
        spec2=species[parts[0]]
        spec2str=parts[0]

        spec1,spec2=min(spec1,spec2),max(spec1,spec2)
        CMid=(spec1,spec2,CMcode)


        newCis=CisMod([spec1str,spec2str],CMcode)
        if CMcode<20:
            cisModules.append(newCis)
        #print line
        
        line=handle.readline() #first malign row
        parts=line.split("\t") 

        while len(line)>0 and parts[2]!="CisModule":
            
            node1Id=(species[parts[0]],parts[2],int(parts[3]))
            node1Line=line

            sPos1=int(parts[3])
            
            line=handle.readline()  #Second malign row
            parts=line.split("\t")

            sPos2=int(parts[3])

            newCis.addSite({spec1str:sPos1,spec2str:sPos2},parts[2],float(parts[5]))
            newCis.sites[-1].addPair(spec1str,spec2str)
                           
            assert(parts[1]=="malign")
            
            line=handle.readline()   #next malign or CisModule
            parts=line.split("\t")
            #print parts[1],



cis3spec=[]
for i in range(len(cisModules)):
    for j in range(i+1,len(cisModules)):
        c3=cisModules[i].intersect(cisModules[j])
        if c3:
            cis3spec.append(c3)

multia=[]
for a in cis3spec:
    for b in cisModules:
        c=a.addPairAlign(b)
        if c:
            multia.append(c)
            print c

#a=cis3spec[0]
#b=cis3spec[16]
#c=a.addPairAlign(cisModules[39])
