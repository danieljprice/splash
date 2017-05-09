#!/bin/bash
#
# short script to copy all splash
# files to a new prefix
#
# ie. splash.defaults, splash.limits, splash.units etc.
# become new.defaults, new.limits, new.units
#
# SPLASH can be invoked to use the new settings files
# using the -p command line option
#
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
   echo 'SPLASH files copy utility -- '
   echo 'copies splash.defaults, splash.limits etc. '
   echo '    to new.defaults, new.limits etc. (use with splash -p new)'
   echo
   echo "Usage $0 newprefix [oldprefix]";
   echo
   echo '(default old prefix is "splash")';
   exit;
else
   new=$1;
   if [ $# -eq 2 ]; then
      old=$2;
   else
      old='splash';
   fi
   for ext in defaults limits units titles anim legend columns filenames; do
       if [ -e $old.$ext ]; then
          cp $old.$ext $new.$ext;
          echo "$old.$ext -> $new.$ext";
       else
          echo "$old.$ext does not exist";
       fi
   done
fi
