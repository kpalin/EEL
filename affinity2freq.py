
import sys

a=sys.stdin.readlines()
b=map(lambda x:map(int,x.split()),a)

tot=[0.0]*len(b[0])

for c in b:
    tot=map(lambda x,y:x+y,c,tot)

suurin=max(tot)
suurin=reduce(lambda x,y:x+y,tot)
c=map(lambda x:map(lambda y,t:round(suurin*y/t),x,tot),b)
d=map(lambda y:"".join(map(lambda x:"%5.0f"%x,y))+"\n",c)


for i in d:
    print i,
