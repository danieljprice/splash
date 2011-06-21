#!/bin/bash
#
# Script supplied by Ben Ayliffe
# for making movies that are playable
# on both Mac and Linux without the need
# for special software
#
# Requires the ffmpeg and mencoder utilities
#
if [ "$1" == "-h" ]
then
    echo '-------------------------------------------------------------------------------'
    echo "Expected input form : <image_files i.e.pgplot> <video name> <quality l/m/h/v/c>"
    echo 'Video quality: [L]ow/[M]ed/[H]igh/[V]eryHigh/[C]ustom'
    echo '<image_files> works as though followed by a wildcard: e.g. pgplot -> pgplot*'
    echo ''
    echo 'Example: newmov pgplot myvideo v'
    echo ''
    echo "If using custom quality option 'c', follow with a bitrate in Kb, 'c 2000'"
    echo '-------------------------------------------------------------------------------'
    exit
fi

vname=`echo $2.avi`
oname=`echo $2.mp4`
qual=$3

#echo 'Video quality: [L]ow/[M]ed/[H]igh/[V]eryHigh/[C]ustom'
#read -e qual

if [ "$qual" == "L" -o "$qual" == "l" ]
then
    qval='-b 300k'
    vval='hq'
    opts='-lavcopts vcodec=mpeg4:vbitrate=300'
fi
if [ "$qual" == "M" -o "$qual" == "m" ]
then
    qval='-b 1000k'
    vval='hq'
    opts='-lavcopts vcodec=mpeg4:vbitrate=1000'
fi
if [ "$qual" == "H" -o "$qual" == "h" ]
then
    qval='-b 5000k'
    vval='hq'
    opts='-lavcopts vcodec=mpeg4:vbitrate=5000'
fi
if [ "$qual" == "V" -o "$qual" == "v" ]
then
    qval=''
    vval='lossless_fast'
    opts='-lavcopts vcodec=ffv1'
fi
if [ "$qual" == "C" -o "$qual" == "c" ]
then
    bnum=$4
    qval=`echo -b $bnum'k'`
    vval='hq'
    opts=`echo -lavcopts vcodec=mpeg4:vbitrate=$bnum`
fi


mencoder mf://$1 -ovc lavc $opts -o $vname
ffmpeg -i $vname -vcodec libx264 -vpre $vval $qval $oname

rm $vname
echo 'Video created: ' $oname
