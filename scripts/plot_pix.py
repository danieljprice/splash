#!/usr/bin/env python
#
# Script to plot the contents of a .pix file output
# with splash -o ascii
#
import numpy as np
import matplotlib.pyplot as plt
import sys
import re

def read_header( file ):
    """
      read the x, y and v plot limits from the .pix header lines
    """
    pat = re.compile(r'.*min\s*=\s+(\S+)\s+max\s*=\s+(\S+)')
    fh = open(file, 'r')
    count = 0
    got = 0
    xmin = np.full((3),0.)
    xmax = np.full((3),1.)
    while (count < 10):
       count += 1

       # Get next line from file
       line = fh.readline()

       # if line is empty
       # end of file is reached
       if not line:
          break

       # otherwise match lines like "min = 0.000 max = 1.000"
       if (pat.match(line)):
          [m] = pat.findall(line)
          xmin[got] = m[0]
          xmax[got] = m[1]
          got += 1

    return xmin[2],xmax[2],xmin[1],xmax[1],xmin[0],xmax[0]

def read_pix( file ):
    """
      read the floating point pixel values
    """
    array = np.loadtxt(file)
    return array

def plot_pix( file ):
    """
      plot the pixel map with appropriate limits
    """
    xmin,xmax,ymin,ymax,vmin,vmax = read_header(file)
    img = read_pix(file)
    print(file,img.shape)
    plt.imshow(img,cmap='RdBu',origin='lower',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax])
    plt.show()

if (len(sys.argv) > 1):
   for file in sys.argv[1:]:
       plot_pix( file )
else:
   print('Usage:',strsys.argv[0],'file.pix')
   sys.exit(2)
