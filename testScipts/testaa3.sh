
source common.sh




## getTFBS
## getTFBSabsolute
## savematch
## showmatch



eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBS aslkjf -showmatch
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBS  -savematch
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBS -savematch /eioikeutta.gff 
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBS 0.9 -showmatch  -savematch haut1.gff -showmatch
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBS -savematch haut2.gff
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBSabsolute 999 -savematch haut3.gff
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBSabsolute 9 -savematch haut4.gff
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBSabsolute  -savematch haut5.gff
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -getTFBSabsolute -999 -savematch haut6.gff

eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG4.bg -getTFBS -savematch haut7.gff
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG4.bg -getTFBSabsolute -savematch haut8.gff


eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -rs ENSRNOG00000006308 SINFRUG00000128278 ENSDARG00000006837 -getTFBSabsolute  -savematch hum_mus.gff



## ALIGN 

eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -rs ENSG00000134323 ENSRNOG00000006308 SINFRUG00000128278 ENSDARG00000006837 -getTFBSabsolute  -align -showalign -savealign -savealignAnchor -savealignGFF
eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -rs ENSRNOG00000006308 SINFRUG00000128278 ENSDARG00000006837 -getTFBSabsolute  -align -showalign -savealign -savealignAnchor -savealignGFF
eelCMD -align hum_mus.gff -showalign -savealign -savealignAnchor -savealignGFF

eelCMD -am FREAC*.pfm -as Nmyc.Enhancers.50kBp.fa -rs ENSRNOG00000006308 SINFRUG00000128278 ENSDARG00000006837 -more -getTFBSabsolute  -align  -more -showalign turha -savealign humMus1.aln -savealignAnchor humMus1.anc -savealignGFF humMus1.gff

eelCMD -align hum_mus.gff -more 100 -showalign -savealign humMus2.aln -savealignAnchor humMus2.anc -savealignGFF humMus2.gff

eelCMD  -more 100 -showalign -savealign humMus2.aln -savealignAnchor humMus2.anc -savealignGFF humMus2.gff

## more
## showalign
## savealign
## savealignAnchor
## savealignGFF




## getTFBS
## getTFBSabsolute
## savematch
## showmatch
