#!/bin/bash
#
# script to make an mpeg4 movie from splash_*.png files
# requires ffmpeg utility
#
# DJP, Feb 2014
#
opts='-r 10 -b 2M -bt 4M'
codec='-vcodec mpeg4'
ffmpeg -i splash_%04d.png $opts $codec movie.mp4
