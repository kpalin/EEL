#!/usr/bin/python
# -*- coding: UTF-8 -*- 


#
# $Log$
#

import re,math

def trunc(val):
    "Round to integer, towards zero"
    rval=int(abs(val))
    if val<0.0:
        rval=-rval

    return rval

class cisMod:
    def __init__(self,id,score):
        self.id=id
        self.score=score

    def setData(self,*dat):

        self.data=dict(dat)
        pass
    


class alnScores:
    "Helper class for analysing proportions for penalty scores"
    def __init__(self,in_str):
        lines=in_str.split("\n")

        assert(lines[0][:3]=="###")
        self.parseHead(lines[0])
        del(lines[0])

        lines=[x.split("\t") for x in lines]
        lines=[x for x in lines if len(x)>7]

        self.parseLines(lines)

        self.splitModules()

    def splitModules(self):
        "Split scores for all of the initialized cis modules"

        self.splits={}
        for mID,mDat in self.cisMods.items():
            #print mID,
            self.splits[mID]=self.splitScores(mDat[1])
        

    def splitScores(self,modData):
        "Split scores for one module"
        aligned=apply(zip,modData.values())


        splitted=[self.splitScore(aligned[0],None)]
        for i in range(1,len(aligned)):
            splitted.append(self.splitScore(aligned[i],aligned[i-1]))

        return splitted


    def squaremod(self,val):
        f=abs(val)
        f-=2*math.pi*trunc(f/(2*math.pi))

        return f*f

    def anglepenalty(self,d,D):
        theta=(d-D)*2.0*math.pi/self.nuc_per_rot

        ret=self.squaremod(theta)/(d+D)

        return ret
            

    def splitScore(self,cur,prev):
        "Split Score for one pair of sites"


        lmbda_p=self.lmbda*(cur[0][3]+cur[1][3])
        xi_p,nu_p,mu_p=0,0,0
        d,D=-1,-1

        try:
            d=cur[0][0]-prev[0][1]-1
            D=cur[1][0]-prev[1][1]-1
            try:
                assert(d>=0)
                assert(D>=0)
            except AssertionError:
                print d
                print D
                raise

            mu_p=self.mu*(d+D)/2.0
            try:
                xi_p=self.xi*self.anglepenalty(d,D)
                nu_p=self.nu*(d-D)**2/float(d+D)
            except ZeroDivisionError:
                pass
        except TypeError:
            pass
        
        score_delta=lmbda_p-xi_p-nu_p-mu_p

        if abs(score_delta-cur[0][2])>0.1:
            print score_delta-cur[0][2],d-D,(lmbda_p,xi_p,nu_p,mu_p)

        return (lmbda_p,xi_p,nu_p,mu_p)





    def parseLines2(self,lines):

        # Get modules
        l=[(x[-1].strip(";"),[float(x[5]),None]) for x in lines if x[2]=="CisModule"]
        self.cisMods=dict(l)

        # Get sequences
        self.seqs={}.fromkeys([x[0] for x in lines]).keys()


        # Get sites for each module
        for mod in self.cisMods.keys():

            # Site data:  (Begin,End,ScoreDelta,MatrixScore)
            
            cmod=dict([(seq, [(int(x[3]),int(x[4]),float(x[5]),float(x[9][3:])) for x in lines if x[2]!="CisModule" and x[0]==seq and x[8]==mod]) for seq in self.seqs])
            self.cisMods[mod][1]=cmod
            
            #print mod,len(cmod),[sum([y[2] for y in x]) for x in cmod.values()]
            
    def parseLines(self,lines):

        # Get modules
        self.cisMods=[cisMod(x[-1].strip(";"),[float(x[5]),None]) for x in lines if x[2]=="CisModule"]

        # Get sequences
        self.seqs={}.fromkeys([x[0] for x in lines]).keys()


        # Get sites for each module
        for mod in self.cisMods:

            # Site data:  (Begin,End,ScoreDelta,MatrixScore)
            
            cmod=apply(mod.setData,[(seq, [(int(x[3]),int(x[4]),float(x[5]),float(x[9][3:])) for x in lines if x[2]!="CisModule" and x[0]==seq and x[8]==mod]) for seq in self.seqs])
            
            #print mod,len(cmod),[sum([y[2] for y in x]) for x in cmod.values()]
            



    def parseHead(self,headLine):
        headRegex=r"### lambda=(?P<lambda>\d+.\d+) mu=(?P<mu>\d+.\d+) nu=(?P<nu>\d+.\d+) xi=(?P<xi>\d+.\d+) Nucleotides per rotation=(?P<nuc_per_rot>\d+.\d+)"
        hr=re.compile(headRegex)
        res=hr.search(headLine)
        assert(res)

        self.lmbda=float(res.group("lambda"))
        self.mu=float(res.group("mu"))
        self.nu=float(res.group("nu"))
        self.xi=float(res.group("xi"))
        self.nuc_per_rot=float(res.group("nuc_per_rot"))

if __name__=="__main__":
    import sys
    a=alnScores(open(sys.argv[1]).read())
    
    

        
