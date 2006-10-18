#
# $Log$
# Revision 1.1.2.2  2006/01/12 10:21:52  kpalin
# Some testing alterations.
#
# Revision 1.1.2.1  2005/03/22 12:21:48  kpalin
# Testing tools.
#
#

BASEDIR=/home/kpalin/tyot/comparative/distEEL/test
PYTHONPATH=$BASEDIR/lib/python2.3/site-packages/:$PYTHONPATH:

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
