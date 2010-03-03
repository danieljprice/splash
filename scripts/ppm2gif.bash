#!/bin/bash
# @(#) converts all ppms to gifs (using Netpbm tools)

for x in splash_*.ppm;
do if [ -e ${x/.ppm/.gif} ]; then
   echo ${x/.ppm/.gif} already exists; 
   else echo creating ${x/.ppm/.gif}; ppmquant 256 $x | ppmtogif > ${x/.ppm/.gif};
   fi;
done;
