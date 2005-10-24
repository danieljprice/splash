#!/bin/bash
# @(#) renames pgplot filenames so they are listed in the correct order
# @(#) (replaces _ with _0 in single then double digit filenames)
# @(#) NB: only works for < 10000 files

for x in pgplot.png_?;
do if test $x != 'pgplot.png_?'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

for x in pgplot.png_??;
do if test $x != 'pgplot.png_??'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

for x in pgplot.png_???;
do if test $x != 'pgplot.png_???'; then echo $x ${x/.png_/_0}.png; mv $x ${x/.png_/_0}.png; fi;
done;
