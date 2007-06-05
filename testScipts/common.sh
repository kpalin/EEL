#
# $Log$
# Revision 1.2  2006/10/18 06:48:13  kpalin
# Added testing scripts.
#
# Revision 1.1.2.2  2006/01/12 10:21:52  kpalin
# Some testing alterations.
#
# Revision 1.1.2.1  2005/03/22 12:21:48  kpalin
# Testing tools.
#
#

BASEDIR=$PWD/../subzero/

echo "BASEDIR=$BASEDIR"

PYTHONPATH=`find $BASEDIR -name site-packages -type d -printf "%p:"`$PYTHONPATH:

EEL=$BASEDIR/bin/eel

OUTFILE=out.tmp
CMPFILE=cmp.tmp

ECHO="echo -e "

eelCMD () 
{
echo "######################"
echo "Running command:"
echo $EEL "$*"
echo -e "######################"
$EEL $*
echo -e "######################\n"
}
