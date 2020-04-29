#!/bin/bash
#
# A basic script to retrieve and install packages to a local directory
# (e.g. in a users home space)
#
# We assume packages can be compiled in the "standard" way using
# "configure" and "make"
#
# A better alternative is to use your inbuilt package manager to install things
#  e.g.
#   Debian/Ubuntu:
#      sudo apt-get install package_name
#   Fedora/Red Hat/CentOS:
#      sudo yum install package_name
#   OpenSUSE:
#      zypper install package_name
#   MacPorts:
#      sudo port install package_name
#   Homebrew:
#      brew install package_name
#
# Written by Daniel Price, April 2020
# Contact: daniel.price@monash.edu
#
xzdist=xz-5.2.1.tar.gz;
xzurl="http://tukaani.org/xz/";
if [ $# -le 1 ]; then
   echo "Usage: $0 http://blah.org/release/blah.tar.gz install_dir";
   exit 1;
fi
disturl=$1;
installprefix=$2;
distfile=$(basename $disturl);
pkg_name=$(basename $distfile .tar.gz)
pkg_name=$(basename $pkg_name .tar.xz)
extension=${distfile/$pkg_name/}
pkg_dir=${distfile/$extension/};
#
#--Check that the install dir is present.
#  This is not strictly necessary, but it means we install cairo and
#  pixman to the same location as the giza libraries and linking of
#  giza with cairo will work automatically.
#
check_install_dir_exists()
{
  if [ ! -d $installprefix ]; then
     echo;
     echo " ERROR: installation directory $installprefix does not exist "
     echo;
     return 1;
  fi
}
#
#--if not already downloaded, retrieve the tarball using wget
#
download_dist_file()
{
  if [ ! -f $distfile ]; then
     echo "$distfile not downloaded";
     if !(type -p wget); then
        echo "ERROR: $0 requires the \"wget\" command, which is not present on";
        echo "your system. Instead, you will need to download the following file by hand:"; echo
        echo "$disturl";
        echo; echo "To proceed, download this files, place them in the current directory and try again"
        return 1;
     else
        wget $disturl;
     fi
  fi
  if [ ! -f $distfile ]; then
     echo; echo "ERROR: $distfile download failed. Please try again"; echo;
     return 1;
  else
     echo "$distfile found in current dir";
     return 0;
  fi
}
#
#--unpack the distribution file with tar or unxz, depending on compression
#
unpack_dist_file()
{
   echo ":: unpacking $distfile to $installprefix";
   if [ $extension==".tar.xz" ]; then
      tar -Jxf $distfile;
   else
      tar xfz $distfile;
   fi
   if [ ! -d $pkg_dir ]; then
      echo; echo "ERROR: failed to unpack (no directory $pkg_dir)"; echo;
      return 1;
      if [ $extension==".tar.xz" ]; then
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
         make || ( echo; echo "ERROR during xzutils build"; echo; return $? );
         make install || ( echo; echo "ERROR installing xzutils into $xzinstalldir"; echo; return $? );
         cd ..;
      #
      #--now unpack using xz utils
      #
         ${xzinstalldir}/bin/unxz $distfile;
         tar xf ${distfile/.xz/};
         if [ ! -d $pkg_dir ]; then
            echo; echo "ERROR: failed to unpack even with xz downloaded (no directory $pkg_dir)"; echo;
            return 1;
         fi
      fi
   fi
}
#
#--install the package to the local file system
#
install_package()
{
   echo ":: installing $pkg_name"
   cd $pkg_dir;
   ./configure --prefix=$installprefix > /dev/null || ( echo; echo "ERROR during config"; echo; return $? );
   make > /dev/null || ( echo; echo "ERROR during build"; echo; return $? );
   make install > /dev/null || ( echo; echo "ERROR installing into $installdir"; echo; return $? );
   cd ..;
}
check_and_finish()
{
   echo ":: $pkg_name installation successful"; echo;
   echo "type \"make\" to compile your main program"; echo;
   echo "You should also add the following line to your .bashrc or equivalent:"; echo;
   if [[ `uname` =~ Darwin ]]; then
      echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$installprefix/lib";
   else
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installprefix/lib";
   fi
   echo;
}
check_install_dir_exists; err=$?;
if [ $err -gt 0 ]; then exit $err; fi
download_dist_file; err=$?;
if [ $err -gt 0 ]; then exit $err; fi
unpack_dist_file; err=$?;
if [ $err -gt 0 ]; then exit $err; fi
install_package; err=$?;
if [ $err -gt 0 ]; then exit $err; fi
check_and_finish; err=$?;
if [ $err -gt 0 ]; then exit $err; fi
