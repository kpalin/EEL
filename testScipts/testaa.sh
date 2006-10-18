
source common.sh

eelCMD  -?
eelCMD -h
eelCMD -help

eelCMD -about -reset

eelCMD -rm 1 -pm -pmw -addMatrix -am eelCMD -addMatrix FREAC* -pm -printMatrices -pmw -printMatrixWeights -removeMatrix 1 -pm -rm -pm -rm 1 -pm -resetMatrices -resm -pm

eelCMD -ps -addSequence -as eelCMD -as Nmyc.Enhancers.50kBp.fa -ps -printSeqNames -rs -rs SINFRUG00000128278 -ps -ress -ps

eelCMD -as invalid.fa -ps

eelCMD -ass Nmyc.Enhancers.50kBp.fa -ps
eelCMD -addSingleSequence invalid.fa -ps

eelCMD -cd

eelCMD -cd .. -cd -dir

eelCMD -q
eelCMD -quit


eelCMD -no-gui<<EOF
quit
EOF


