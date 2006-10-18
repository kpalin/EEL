
source common.sh


eelCMD "-setMarkovBG" 

eelCMD -setMarkovBG $EEL

eelCMD -setMarkovBG ENSRNOG00000006308

eelCMD -as  Nmyc.Enhancers.50kBp.fa -ps -setMarkovBG SINFRUG00000128278 ENSRNOG00000006308

eelCMD -as  Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG00000128278 -saveMarkovBackground SINFRUGdefault.bg
eelCMD -as  Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG00000128278 3 -saveMarkovBackground SINFRUG3.bg
eelCMD -as  Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG00000128278 4 -saveMarkovBackground SINFRUG4.bg
eelCMD -as  Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG00000128278 5 -saveMarkovBackground SINFRUG.bg
eelCMD -as  Nmyc.Enhancers.50kBp.fa -setMarkovBG SINFRUG00000128278 6 -saveMarkovBackground SINFRUG6.bg

eelCMD -setMarkovBG SINFRUG3.bg
eelCMD -setMarkovBG SINFRUG4.bg
eelCMD -setMarkovBG SINFRUG5.bg
eelCMD -setMarkovBG SINFRUG6.bg

eelCMD -am FREAC* -pmw
eelCMD -setMarkovBG SINFRUG.bg -am FREAC* -pmw
eelCMD -am FREAC* -setMarkovBG SINFRUG.bg  -pmw
eelCMD -am FREAC* -setMarkovBG SINFRUG.bg  -pmw -setBGfreq a b c d e f -pmw
eelCMD -am FREAC* -setMarkovBG SINFRUG.bg  -pmw -setBGfreq 0.25 0.25 0.25 0.25 -pmw

eelCMD -cd .. -am */FREAC* -setMarkovBG */SINFRUG.bg  -pmw -setBGfreq 0.25 0.25 0.25 0.25 -pmw

eelCMD -setpseudocount 0.01 -am FREAC* -setMarkovBG SINFRUG.bg  -pmw -setBGfreq 0.25 0.25 0.25 0.25 -pmw
eelCMD  -am FREAC* -setpseudocount 0.01 -setMarkovBG SINFRUG.bg  -pmw -setBGfreq 0.25 0.25 0.25 0.25 -pmw
eelCMD -am FREAC* -setMarkovBG SINFRUG.bg  -setpseudocount 0.01 -pmw -setBGfreq 0.25 0.25 0.25 0.25 -pmw
eelCMD  -am FREAC* -setMarkovBG SINFRUG.bg  -pmw -setBGfreq 0.25 0.25 0.25 0.25 -setpseudocount 0.01 -pmw



## ALIGN 
## more
## showalign
## savealign
## savealignAnchor
## savealignGFF

## getTFBS
## getTFBSabsolute
## savematch
## showmatch
