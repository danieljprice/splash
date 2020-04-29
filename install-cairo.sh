#!/bin/bash
#
# Script for splash 2.x that retrieves and installs
# both cairo and pixman
#
# (these are the only dependencies for the giza backend,
#  are often already present as system libraries but
#  may need to be installed by the user if not)
#
# An alternative is to use your inbuilt package manager to install cairo
#  e.g.
#   Debian/Ubuntu:
#      sudo apt-get install libcairo2-dev
#   Fedora/Red Hat/CentOS:
#      sudo yum install cairo-devel
#   OpenSUSE:
#      zypper install cairo-devel
#   MacPorts:
#      sudo port install cairo
#
installprefix=$PWD/giza;
./install-pkg.sh "http://cairographics.org/releases/pixman-0.40.0.tar.gz" $installprefix
./install-pkg.sh "http://cairographics.org/releases/cairo-1.16.0.tar.xz" $installprefix
echo "type \"make\" to compile SPLASH"; echo;
