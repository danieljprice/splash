#!/bin/bash
#
# A basic script to retrieve and install packages to a local directory
# (e.g. in a users home space)
#
# We assume packages can be compiled in the "standard" way using
# "configure" and "make"
#
# Usage: install-pkg.sh <url> <prefix> [configure_extra...]
#
# Written by Daniel Price, April 2020
# Contact: daniel.price@monash.edu
#
xzdist=xz-5.2.1.tar.gz;
xzurl="http://tukaani.org/xz/";
if [ $# -lt 2 ]; then
   echo "Usage: $0 <url> <install_dir> [configure_extra...]";
   exit 1;
fi
disturl=$1;
installprefix=$2;
shift 2
configure_extra="$*"
distfile=$(basename $disturl);
pkg_name=$(basename $distfile .tar.gz)
pkg_name=$(basename $pkg_name .tar.xz)
extension=${distfile/$pkg_name/}
pkg_dir=${distfile/$extension/};
#
#--parallel make: inherit MAKEFLAGS from the environment when set;
# otherwise pass an explicit -j flag to make.
#
if [ -z "${MAKEFLAGS:-}" ]; then
   nproc_val=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
   make_j="-j${nproc_val}"
else
   make_j=""
fi
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
#--if not already downloaded, retrieve the tarball
#
download_dist_file()
{
  if [ ! -f $distfile ]; then
     echo "$distfile not downloaded";
     if type -p wget > /dev/null 2>&1; then
        wget $disturl;
     elif type -p curl > /dev/null 2>&1; then
        curl -LO $disturl;
     else
        echo "ERROR: $0 requires wget or curl, which is not present on your system.";
        echo "Please download the following file by hand:"; echo
        echo "$disturl";
        return 1;
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
   if [ "$extension" = ".tar.xz" ]; then
      tar -Jxf "$distfile";
   else
      tar xfz "$distfile";
   fi
   if [ ! -d "$pkg_dir" ]; then
      if [ "$extension" != ".tar.xz" ]; then
         echo; echo "ERROR: failed to unpack (no directory $pkg_dir)"; echo;
         return 1;
      else
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
         make $make_j || ( echo; echo "ERROR during xzutils build"; echo; return $? );
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
   cd "$pkg_dir" || return 1
   if [ ! -f ./configure ]; then
      echo; echo "ERROR: no ./configure in $pkg_dir (autotools required)"; echo
      cd ..; return 1
   fi
   if [ -n "$CC" ]; then
      export CC
   fi
   if [ -n "$CFLAGS" ]; then
      export CFLAGS
   fi
   if [ -n "$PKG_CONFIG_PATH" ]; then
      export PKG_CONFIG_PATH
   fi
   echo ":: configure $pkg_name (CC=${CC:-default} CFLAGS=${CFLAGS:-default})"
   if ! ./configure --prefix="$installprefix" $configure_extra; then
      echo; echo "ERROR during config for $pkg_name"; echo
      cd ..; return 1
   fi
   if ! make $make_j; then
      echo; echo "ERROR during build for $pkg_name"; echo
      cd ..; return 1
   fi
   if ! make install; then
      echo; echo "ERROR installing $pkg_name into $installprefix"; echo
      cd ..; return 1
   fi
   cd ..
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
