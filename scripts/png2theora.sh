#!/bin/sh
#
# Converts sequence of pgplot.png* files to Ogg Theora movie
# using the png2theora utility
#
# Contributed by Pau Amaro-Seoane (pau@aei.mpg.de)
# May 2011
#
i=1;
temp=$(mktemp -p .);
for file in $(ls pgplot.png* | sort -V); do
    mv "$file" $temp;
    mv $temp $(printf "image_%0.3d.png" $i);
    i=$((i + 1));
done

png2theora --video-quality 3 --framerate-numerator 5 \
          --framerate-denominator 1 --keyframe-freq 12 \
          --chroma-444 \
                               image_%03d.png \
          --output output.ogg
