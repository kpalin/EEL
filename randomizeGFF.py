import random,sys


fname=sys.argv[1]

lines=map(lambda x:x.split("\t"),open(fname).readlines())

names=map(lambda x:(x[2],int(x[4])-int(x[3])),lines)

random.shuffle(names)


m={}
goodLines=[]
for (name,dim) in names:
    line=lines.pop()
    c=0
    while m.has_key(line[3]) and m[line[3]].has_key(name):
        c=c+1
        lines=[line]+lines
        prevLine=line
        line=lines.pop()
        #print "looping with",line
        if c>10:
            random.shuffle(lines)
    if not m.has_key(line[3]):
        m[line[3]]={}
    m[line[3]][name]=1
    line[2]=name
    line[4]=str(int(line[3])+dim)
    #print line
    goodLines.append(line)


#print goodLines
for x in goodLines: print "\t".join(x),

sys.exit()

for line,(name,dim) in zip(lines,names):
    line[2]=name
    line[4]=str(int(line[3])+dim)

#lines=map(lambda line:"%s\n"%("\t".join(line)),lines)
lines.sort(lambda x,y:cmp((x[2],int(x[3])),(y[2],int(y[3]))))
           
for i in range(len(lines)):
    if i>0 and lines[i][2]==lines[i-1][2] and lines[i][3]==lines[i-1][3] and lines[i][6]==lines[i-1][6]:
        #print "Skipping",i,"\t".join(lines[i])
        continue
    print "\t".join(lines[i]),
    
