from matrix import *
from mmap import *

#a=open("/home/kpalin/tyot/comparative/humanGenome/chr21.fa")
a=open("/home/kpalin/tyot/comparative/mabs/random.fa")
#a=open("../lyhyt.fa")
a.seek(0,2)
l=a.tell()
a.seek(0)

#bgm=mmap(a.fileno(),l,prot=PROT_READ,access=ACCESS_READ)
bgm=mmap(a.fileno(),l,prot=PROT_READ)

print "Computing background"

counts={}
#for i in bgm:
#    if not counts.has_key(i):
#        counts[i]=0
#    counts[i]+=1

bg=BackGround(a)

mattup=bg.giveGramVector()
print mattup

bg2=BackGround(mattup)
#print bg2.giveGramVector()

print "Done counting"
grams=bg.giveGramCounts().items()
print "Got counts"
grams.sort()



#for k,v in grams:
#    print "%s %d\t %g"%(k,v,bg.stringProb(k))
