#!/usr/bin/python
# -*- coding: UTF-8 -*- 


#
# $Log$
# Revision 1.3  2005/01/03 11:25:00  kpalin
# Works flawlessly.
#
# Revision 1.2  2004/12/30 11:31:52  kpalin
# Classified.
#
# Revision 1.1  2004/12/30 11:11:54  kpalin
# Working version.
#
#

import re,math
from cStringIO import StringIO


def trunc(val):
    "Round to integer, towards zero"
    rval=int(abs(val))
    if val<0.0:
        rval=-rval

    return rval


class alignedSite:
    def __init__(self,settings,alnThis,alnPrev=None):
        self.weight1,self.weight2=alnThis[0][3],alnThis[1][3]
        self.d,self.D=-1,-1
        self.report_score=alnThis[0][2]
        self.param=settings

        try:
            self.d=alnThis[0][0]-alnPrev[0][1]-1
            self.D=alnThis[1][0]-alnPrev[1][1]-1
            try:
                assert(self.d>=0 and self.D>=0)
            except AssertionError:
                print self.d,self.D
                raise
        except TypeError:
            pass
        
        self.splitScore(alnThis,alnPrev)

        
    def squaremod(self,val):
        f=abs(val)
        f-=2*math.pi*trunc(f/(2*math.pi))

        return f*f

    def anglepenalty(self,d,D):
        theta=(d-D)*2.0*math.pi/self.param.nuc_per_rot

        ret=self.squaremod(theta)/(d+D)

        return ret
            

    def splitScore(self,cur,prev):
        "Split Score for one pair of sites"


        self.lmbda_p=self.param.lmbda*(self.weight1+self.weight2)
        self.xi_p,self.nu_p,self.mu_p=0,0,0

        if(self.d>0 or self.D>0):
            self.mu_p=self.param.mu*(self.d+self.D)/2.0
            try:
                self.xi_p=self.param.xi*self.anglepenalty(self.d,self.D)
                self.nu_p=self.param.nu*(self.d-self.D)**2/float(self.d+self.D)
            except ZeroDivisionError:
                pass
        
        self.score_delta=self.lmbda_p-self.xi_p-self.nu_p-self.mu_p

        if abs(self.score_delta-self.report_score)>0.1:
            print self.score_delta-self.report_score,self.d,self.D,(self.lmbda_p,self.xi_p,self.nu_p,self.mu_p)



    def __str__(self):
        s="%6.2f = %5.2f - %5.2f(%3.2f) - %5.2f(%d) - %5.2f = weight - dist - deltaDist - rot"%(self.score_delta,self.lmbda_p,self.mu_p,(self.d+self.D)/2.0,self.nu_p,abs(self.d-self.D),self.xi_p)
        return s
    
def mean(lst):
    return sum(lst)*1.0/len(lst)

def varAndMean(lst,meanVal=None):
    if not meanVal:
        meanVal=mean(lst)
    varVal=mean([(x-meanVal)**2 for x in lst])
    return (varVal,meanVal)

def var(lst,meanVal=None):
    varVal,meanVal=varAndMean(lst,meanVal)
    return varVal

        

class cisMod:
    def __init__(self,id,score,param):
        self.id=id
        self.score=score
        self.param=param

    def setData(self,*dat):
        self.data=dict(dat)
        pass

    def termsStat(self):
        self.lmbda_var,self.lmbda_mean=varAndMean([x.lmbda_p for x in self.splitted])
        self.mu_var,self.mu_mean=varAndMean([x.mu_p for x in self.splitted])
        self.nu_var,self.nu_mean=varAndMean([x.nu_p for x in self.splitted])
        self.xi_var,self.xi_mean=varAndMean([x.xi_p for x in self.splitted])

    def splitScores(self):
        "Split scores for the module"
        aligned=apply(zip,self.data.values())


        self.splitted=[alignedSite(self.param,aligned[0],None)]
        self.splitted.extend([alignedSite(self.param,aligned[i],aligned[i-1]) for i in range(1,len(aligned))])


    def __str__(self):
        if not hasattr(self,"lmbda_mean"):
            self.termsStat()
        a=StringIO()
        a.write("%s: %f(+-%f)-%f(+-%f)-%f(+-%f)-%f(+-%f)\n"%(self.id,self.lmbda_mean,math.sqrt(self.lmbda_var),\
                                                             self.mu_mean,math.sqrt(self.mu_var),\
                                                             self.nu_mean,math.sqrt(self.nu_var),\
                                                             self.xi_mean,math.sqrt(self.xi_var)))
        a.write("\n".join([str(site) for site in self.splitted]))

        return a.getvalue()


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
        for mod in self.cisMods:
            #print mID,
            mod.splitScores()
        


    def __str__(self):
        a=StringIO()
        a.write("### lambda=%f mu=%f nu=%f xi=%f Nucleotides per rot=%f\n"%(self.lmbda,self.mu,self.nu,self.xi,self.nuc_per_rot))
        a.write("\n".join([str(mod) for mod in self.cisMods]))
        return a.getvalue()


    def parseLines(self,lines):

        # Get modules
        l=[(x[-1].strip(";"),float(x[5]),self) for x in lines if x[2]=="CisModule"]
        l={}.fromkeys(l).keys()
        l.sort(lambda x,y:cmp(y[1],x[1]))
        self.cisMods=[apply(cisMod,x) for x in l]

        # Get sequences
        self.seqs={}.fromkeys([x[0] for x in lines]).keys()


        # Get sites for each module
        for mod in self.cisMods:

            # Site data:  (Begin,End,ScoreDelta,MatrixScore)
            cmod=apply(mod.setData,[(seq, [(int(x[3]),int(x[4]),float(x[5]),float(x[9][3:])) for x in lines if x[2]!="CisModule" and x[0]==seq and x[8]==mod.id]) for seq in self.seqs])
            
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
    for i in a.cisMods[:10]:print i
    print "\n\n\nHUONOJA\n\n\n"
    for i in a.cisMods[-10:]:print i
    

        
