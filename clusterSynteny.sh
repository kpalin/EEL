#!/bin/sh

FSHOME="/fs/home/kpalin/"

#
# $Log$
# Revision 1.2  2004/02/23 12:22:57  kpalin
# Updates for per gene orthologous runs.
#
# Revision 1.1  2004/02/20 11:05:02  kpalin
# More or less the script used to run the syntenic region alignments.
# Added timing
#
#


#

if [ -e $FSHOME ]; then
    rm -rf /tmp/kpalin/mabsrun/
    mkdir /tmp/kpalin/mabsrun/
    cd /tmp/kpalin/mabsrun/

    tar -zxvf $FSHOME/tyot/mabs/mabs-1_beta12.tar.gz
fi


cd mabs-1_beta12
echo $PWD

make

SYNTENYDIR="$FSHOME/tyot/mabs/synteny/perGene/"


RUNID="mabs$$"

FASTAFILE="og_ENSG00000164161_ENSMUSG00000031718.fasta"

SEQFILE="$SYNTENYDIR$FASTAFILE"
echo "Testing alignment with $SEQFILE"
#/usr/bin/time /bin/nice -5 mabs -as $SEQFILE -am TFv1/*.pfm TFv1/*/*.pfm -setMarkovBG human.ChrI.O4.bg -getTFBSabsolute 9 -savematch $FASTAFILE.gff.gz -align $FASTAFILE.gff.gz 200 2.0 200.0 0.12 200.0 . -savealignGFF $FASTAFILE.align.gff 2>&1  >$FASTAFILE.mabs.log 
echo done $FASTAFILE

mv $FASTAFILE.* $SYNTENYDIR/output/

#FASTAFILE=`python2.2 $FSHOME/tyot/mabs/ogRun/nextFasta.py $RUNID`
FASTAFILE=`python2.2 nextFasta.py $RUNID`



# Limit Reserved Set Size to about 1.3GB
ulimit -m 1372160
# Limit Virtual Memory Size to about 1.66GB
ulimit -v $[1024*1024*5/3]

echo $FASTAFILE

while [ "$FASTAFILE" != "xAllDone" ] ; do
    SEQFILE="$SYNTENYDIR$FASTAFILE"
    echo "Starting" $SEQFILE
    /usr/bin/time /bin/nice -5 mabs -as $SEQFILE -am TFv1/*.pfm TFv1/*/*.pfm \
	-setMarkovBG human.ChrI.O4.bg -getTFBSabsolute 9 \
	-align . 200 2.0 200.0 0.12 200.0 . \
	-savealignGFF $FASTAFILE.align.gff 2>&1  >$FASTAFILE.mabs.log 
#    /usr/bin/time /bin/nice -5 mabs -as $SEQFILE -am TFv1/*.pfm TFv1/*/*.pfm \
#	-setMarkovBG human.ChrI.O4.bg -getTFBSabsolute 9 \
#	-savematch $FASTAFILE.gff.gz \
#	-align $FASTAFILE.gff.gz 200 2.0 200.0 0.12 200.0 . \
#	-savealignGFF $FASTAFILE.align.gff 2>&1  >$FASTAFILE.mabs.log 
    gzip $FASTAFILE.mabs.log
    if [ ! $? == 0 ]; then
	echo "GZIP failed"
	exit
    fi

    gzip $FASTAFILE.align.gff
    if [ ! $? == 0 ]; then
	echo "GZIP failed"
	exit
    fi
    echo done $FASTAFILE
    mv $FASTAFILE.* $SYNTENYDIR/output/
    if [ -e $SYNTENYDIR/output/$FASTAFILE.align.gff.gz ]; then
	FASTAFILE=`python2.2 nextFasta.py $RUNID`
    else
	echo "Can't find the result file."
	FASTAFILE="xAllDone"
    fi
    #FASTAFILE=`python2.2 $FSHOME/tyot/mabs/ogRun/nextFasta.py $RUNID`
done
#nice -5 python2.2 ParameterOptimize.py  $*  2> ../ParamOpt.`hostname`.`date +%d%m%y%k%m`.err >/dev/null </dev/null &

