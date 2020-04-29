#!/bin/bash
#
# Script for splash 2.x that retrieves and installs
# cfitsio
#
# An alternative is to use your inbuilt package manager
#
installprefix=$PWD/giza;
./install-pkg.sh "https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz" $installprefix
