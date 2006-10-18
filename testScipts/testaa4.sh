
source common.sh

# suboptimalsDownTo 


eelCMD -align   hum_mus.gff 0 3 1.0 0.05 2.0 -suboptimalsDownTo 150 \
    -showalign -savealign humMus2.aln -savealignAnchor humMus2.anc \
    -savealignGFF humMus2.gff

