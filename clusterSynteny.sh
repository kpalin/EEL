#!/bin/sh

FSHOME="/fs/home/kpalin/"

#
# $Log$
#


#rm -rf /tmp/kpalin/mabsrun/
mkdir /tmp/kpalin/mabsrun/
cd /tmp/kpalin/mabsrun/


tar -zxvf $FSHOME/tyot/mabs/mabs-1_beta12.tar.gz


cd mabs-1_beta12
echo $PWD

make

SYNTENYDIR="$FSHOME/tyot/mabs/synteny/"


RUNID="mabs$$"

FASTAFILE="fasta_hum_mus_1031"

SEQFILE="$SYNTENYDIR$FASTAFILE"
echo "Testing alignment with $SEQFILE"
time nice -5 mabs -as $SEQFILE -am TFv1/*.pfm TFv1/*/*.pfm -setMarkovBG human.ChrI.O4.bg -getTFBSabsolute 9 -savematch $FASTAFILE.gff.gz -align $FASTAFILE.gff.gz 200 2.0 200.0 0.12 200.0 . -savealignGFF $FASTAFILE.align.gff 2>&1  >$FASTAFILE.mabs.log 
echo done $FASTAFILE

mv $FASTAFILE.* $SYNTENYDIR/output/

FASTAFILE=`python2.2 $FSHOME/tyot/mabs/nextFasta.py $RUNID`



# Limit Reserved Set Size to about 1.3GB
ulimit -m 1372160
# Limit Virtual Memory Size to about 1.66GB
ulimit -v $[1024*1024*5/3]

while [ "$FASTAFILE" != "xAllDone" ] ; do
    SEQFILE="$SYNTENYDIR$FASTAFILE"
    echo $SEQFILE
    time nice -5 mabs -as $SEQFILE -am TFv1/*.pfm TFv1/*/*.pfm -setMarkovBG human.ChrI.O4.bg -getTFBSabsolute 9 -savematch $FASTAFILE.gff.gz -align $FASTAFILE.gff.gz 200 2.0 200.0 0.12 200.0 . -savealignGFF $FASTAFILE.align.gff 2>&1  >$FASTAFILE.mabs.log 
    gzip $FASTAFILE.mabs.log
    gzip $FASTAFILE.align.gff
    echo done $FASTAFILE
    mv $FASTAFILE.* $SYNTENYDIR/output/
    FASTAFILE=`python2.2 $FSHOME/tyot/mabs/nextFasta.py $RUNID`
done
#nice -5 python2.2 ParameterOptimize.py  $*  2> ../ParamOpt.`hostname`.`date +%d%m%y%k%m`.err >/dev/null </dev/null &

