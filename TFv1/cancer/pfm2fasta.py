import sys


fname=sys.argv[1]

M=map(lambda x:map(int,x.split()),open(fname).readlines())

n=len(M[0])
sigma=len(M)

lines=0
for i in range(sigma):
    lines+=M[i][0]

chrs=[]
for i in range(lines):
    chrs.append([])

alphabet=["A","C","G","T"]
for i in range(n):
    l=0
    for j in range(sigma):
        for k in range(M[j][i]):
            try:
                chrs[l].append(alphabet[j])
            except Exception:
                print "l=%d j=%d lines=%d"%(l,j,lines)
                raise
            l+=1
chrs=map("".join,chrs)

print "> %s\n%s"%(fname,"\n".join( chrs))
