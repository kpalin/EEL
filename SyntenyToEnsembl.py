
import re,fileinput,sys

# $Log$

Species={"HUMAN":"http://www.ensembl.org/Homo_sapiens/contigview?chr=%d&vc_start=%d&vc_end=%d",
         "MOUSE":"http://www.ensembl.org/Mus_musculus/contigview?chr=%d&vc_start=%d&vc_end=%d"}

m=re.compile(r"(%s)\|(\d+)[.](\d+)-(\d+) "%("|".join(Species.keys())))

def SyntenyToEnsembl(string):
    mat=m.search(string)
    if not mat:
        return None
    else:
        return Species[mat.group(1)]%(int(mat.group(2)),int(mat.group(3)),int(mat.group(4)))
    
def SyntenyLen(string):
    mat=m.search(string)
    if not mat:
        return None
    else:
        return int(mat.group(4))-int(mat.group(3))

for line in fileinput.input():
    link=SyntenyToEnsembl(line)
    if not link:continue
    print link
    print >> sys.stderr, "%s:%d"%(fileinput.filename(),SyntenyLen(line))
