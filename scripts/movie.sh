#!/bin/bash
#
# script to make an mpeg4 movie from png files
# by default, takes splash_*.png as inputs and movie.mp4 as output
# requires ffmpeg utility
#
# DJP, Feb 2014. edited by KAH, Aug 2017
#
opts='-r 10 -vb 50M -bt 100M -vcodec mpeg4 -vf setpts=4.*PTS'
if [ $# -le 0 ]; then
   ffmpeg -i splash_%04d.png $opts movie.mp4
fi
if [ $# -eq 1 ]; then
   ffmpeg -i $1_%04d.png $opts movie.mp4
fi
if [ $# -eq 2 ]; then
   ffmpeg -i $1_%04d.png $opts $2.mp4
fi
if [ $# -ge 3 ]; then
   echo "usage: $0 infile_prefix outfile_prefix"
   exit
fi
