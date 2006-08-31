#!/bin/sh

EELHOME=/fs-3/e/kpalin/tyot/EEL
export PYTHONPATH=$EELHOME/subzero/lib/python2.4/site-packages/


if [ ! -r  $EELHOME/testScipts/test2d.gff ] ; then
    $EELHOME/subzero/bin/eel -as $EELHOME/Nmyc.Enhancer.50kBp.ENSG00000134323.fa -as $EELHOME/Nmyc.Enhancer.50kBp.ENSMUSG00000037169.fa -am $EELHOME/TFv1/*.pfm -getTFBSabsolute 9 -savematch $EELHOME/testScipts/test2d.gff -align -more 50 -savealignGFF $EELHOME/testScipts/test2d.aln.gff 
fi

if [ ! -r  $EELHOME/testScipts/test.rat.gff ] ; then
    $EELHOME/subzero/bin/eel -as $EELHOME/Nmyc.Enhancer.50kBp.ENSRNOG00000006308.fa -am $EELHOME/TFv1/*.pfm -getTFBSabsolute 9 -savematch $EELHOME/testScipts/test.rat.gff  
fi


if [ ! -r  $EELHOME/testScipts/test2d.lyhyt.aln.gff ] ; then
    $EELHOME/subzero/bin/eel -align $EELHOME/testScipts/test2d.lyhyt.gff -savealignGFF $EELHOME/testScipts/test2d.lyhyt.aln.gff
fi


#(gdb) p p.output()
#(295,197[MZF1-4.pfm+])$6 = void
#(gdb) p p.getValue()
#$7 = -39.774248
#(gdb)   

#FOR GDB:  run ../subzero/bin/eel -multipleAlign test2d.lyhyt.gff
#/scratch/kpalin/bin/valgrind --tool=callgrind \
#/scratch/kpalin/bin/valgrind --tool=memcheck  --suppressions=/fs-3/e/kpalin/tyot/EEL/testScipts/p.supp  \
#    python $EELHOME/subzero/bin/eel \
#	-multipleAlign $EELHOME/testScipts/test3d.lyhyt.gff \
#	-savealignGFF $EELHOME/testScipts/test3d.lyhyt.aln.gff &

#	-savealignGFF $EELHOME/testScipts/testMd.lyhyt.aln.gff
#	-multipleAlign $EELHOME/testScipts/test2d.lyhyt.gff \
#	-savealignGFF $EELHOME/testScipts/testMd.lyhyt.aln.gff
#	-multipleAlign $EELHOME/testScipts/test3d.lyhyt.gff \
#	-savealignGFF $EELHOME/testScipts/test3d.lyhyt.aln.gff
#	-align $EELHOME/testScipts/test2d.lyhyt.gff \
#	-savealignGFF $EELHOME/testScipts/testMd.lyhyt.aln.gff
/scratch/kpalin/bin/valgrind --tool=callgrind --simulate-cache=yes \
    python $EELHOME/subzero/bin/eel \
	-multipleAlign $EELHOME/testScipts/test2d.lyhyt.gff \
	-savealignGFF $EELHOME/testScipts/testMd.lyhyt.aln.gff

