#!/bin/bash
# @(#) renames pgplot filenames so they are listed in the correct order
# @(#) NB: only works for < 10000 files
#
# 25/7/2006 Daniel Price, University of Exeter
#           dprice@astro.ex.ac.uk
#
# Usage: fixpgplotnames.bash [noffset]
#        where noffset changes the starting number (default is zero)
#
if [ $# -ne 1 ]; then
   numoffset=0;
else
   numoffset=$1;
   echo Starting file numbers at $numoffset;
   # should check if offset is really a number, otherwise might get garbage
   ##if [ ${numoffset} -ne '[0-9]' ]; then exit "Usage: fixgifs.bash [nstart]"; fi
fi
#
# copy first file(s) (e.g. pgplot.gif) with offset if appropriate
#
for x in pgplot.???;
do
   num=1;
   let "num=num+numoffset";
   lennum=${#num};
   newname=${x/./_$num.}
   if test $lennum -eq 1; then newname=${newname/_/_000}; fi;
   if test $lennum -eq 2; then newname=${newname/_/_00}; fi;
   if test $lennum -eq 3; then newname=${newname/_/_0}; fi;
   echo $x '->' $newname;
   mv $x $newname;
done;
#
# fix all subsequent files (e.g. pgplot.gif_1 pgplot.gif_2 ... pgplot.gif_11 )
#
for x in pgplot.???_*;
do
 num=${x##pgplot*_}; # extract number from the end of the string
 prefix=${x%%_*};    # extract string before the number
#
# add the offset
#
 let "num=num+numoffset";
#
# construct new filename
#
 newname=${prefix/./_$num.};
 lennum=${#num};
#
# add appropriate number of zeros
#
 if test $lennum -eq 1; then newname=${newname/_/_000}; fi;
 if test $lennum -eq 2; then newname=${newname/_/_00}; fi;
 if test $lennum -eq 3; then newname=${newname/_/_0}; fi;
 echo $x '->' $newname; 
 mv $x $newname;
done;
exit;
