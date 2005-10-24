#!/bin/bash
# @(#) renames pgplot filenames so they are listed in the correct order
# @(#) (replaces _ with _0 in single then double digit filenames)
# @(#) NB: only works for < 10000 files

for x in pgplot.gif_?;
do if test $x != 'pgplot.gif_?'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

for x in pgplot.gif_??;
do if test $x != 'pgplot.gif_??'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

for x in pgplot.gif_???;
do if test $x != 'pgplot.gif_???'; then echo $x ${x/.gif_/_0}.gif; mv $x ${x/.gif_/_0}.gif; fi;
done;
