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
cairodist=cairo-1.12.18.tar.xz;
pixmandist=pixman-0.32.6.tar.gz;
xzdist=xz-5.2.1.tar.gz;
installprefix=$PWD/giza;
url="http://cairographics.org/releases";
xzurl="http://tukaani.org/xz/";
#
#--Check that the giza directory is present.
#  This is not strictly necessary, but it means we install cairo and
#  pixman to the same location as the giza libraries and linking of 
#  giza with cairo will work automatically.
#
if [ ! -d $installprefix ]; then
   echo;
   echo " ERROR: directory $installprefix does not exist "; echo;
   echo " $0 should be run from the root-level splash directory";
   echo " with giza already downloaded as a subdirectory of splash";
   echo;
   exit 1;
fi
#
#--if not already downloaded, retrieve the pixman and cairo tarballs using wget
#
if [ ! -f $cairodist ] || [ ! -f $pixmandist ]; then
   echo "cairo and/or pixman not downloaded";
   if !(type -p wget); then
      echo "ERROR: $0 requires the \"wget\" command, which is not present on";
      echo "your system. Instead, you will need to download the following files by hand:"; echo
      echo "$url/$cairodist";
      echo "$url/$pixmandist";
      echo; echo "To proceed, download these files, place them in the current directory and try again"
      exit;
   else
      wget $url/$pixmandist;
      wget $url/$cairodist;
   fi
fi
#
#--proceed with installation
#
if [ ! -f $cairodist ] || [ ! -f $pixmandist ]; then
   echo; echo "ERROR: cairo and/or pixman download failed. Please try again"; echo;
else
   echo "$pixmandist and $cairodist found in current dir";
#
#--unpack the distribution files
#
   echo "unpacking pixman...";
   tar xfz $pixmandist;
   echo "unpacking cairo...";
   tar -Jxf $cairodist;
   pixmandir=${pixmandist/.tar.gz/};
   cairodir=${cairodist/.tar.xz/};
   if [ ! -d $pixmandir ]; then
      echo; echo "ERROR: pixman failed to unpack (no directory $pixmandir)"; echo;
      exit $?;
   fi
   if [ ! -d $cairodir ]; then
      echo; echo "ERROR: cairo failed to unpack (no directory $cairodir)"; echo;
      #
      #--install xzutils if tar -Jxf fails...
      #
      echo "Attempting to download xzutils in order to unpack cairo..."
      wget $xzurl/$xzdist;
      tar xfz $xzdist;
      xzdir=${xzdist/.tar.gz/};
      cd $xzdir;
      xzinstalldir=/tmp/xz-tmp/;
      ./configure --prefix=$xzinstalldir;
      make || ( echo; echo "ERROR during xzutils build"; echo; exit $? );
      make install || ( echo; echo "ERROR installing xzutils into $xzinstalldir"; echo; exit $? );
      cd ..;
      #
      #--now unpack cairo using xz utils
      #
      ${xzinstalldir}/bin/unxz $cairodist;
      tar xf ${cairodist/.xz/};
      if [ ! -d $cairodir ]; then
         echo; echo "ERROR: cairo failed to unpack even with xz downloaded (no directory $cairodir)"; echo;
         exit $?;
      fi
   fi
#
#--install pixman
#
   cd $pixmandir;
   ./configure --prefix=$installprefix || ( echo; echo "ERROR during pixman config"; echo; exit $? );
   make || ( echo; echo "ERROR during pixman build"; echo; exit $? );
   make install || ( echo; echo "ERROR installing pixman into $installdir"; echo; exit $? );
   cd ..;
#
#--install cairo
#
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$installprefix;
   export PKG_CONFIG_PATH=$installprefix/lib/pkgconfig;
   cd $cairodir;
   ./configure --prefix=$installprefix || ( echo; echo "ERROR during cairo config"; echo; exit $? );
   make || ( echo; echo "ERROR during cairo build"; echo; exit $? );
   make install || ( echo; echo "ERROR installing cairo into $installdir"; echo; exit $? );
   cd ..;
#--cleanup
   #rm -rf $pixmandir;
   #rm -rf $cairodir;
#
#--finish
#
   echo; echo "cairo and pixman installation successful";
   echo "type \"make\" to compile SPLASH"; echo;
   echo "You should also add the following line to your .bashrc or equivalent:"; echo;
   if [[ `uname` =~ Darwin ]]; then
      echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$installprefix";   
   else
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installprefix";
   fi
   echo
fi



