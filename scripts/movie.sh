#!/bin/bash
#
# script to make an mpeg4 movie from splash_*.png files
# requires ffmpeg utility
#
# DJP, Feb 2014
#
opts='-r 10 -vb 50M -bt 100M -vcodec mpeg4 -vf setpts=4.*PTS'
ffmpeg -i splash_%04d.png $opts movie.mp4
