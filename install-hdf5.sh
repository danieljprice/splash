#!/bin/bash
#
# Script for SPLASH that retrieves and installs a static HDF5 library
# using install-pkg.sh
#
# Usage: install-hdf5.sh [--prefix DIR] [--disable-shared]
#
# Written by Daniel Price
#
set -e

hdf5_url="https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz"
prefix="${HOME}/local/hdf5"
shared_flag=""

while [ $# -gt 0 ]; do
   case "$1" in
      --prefix)
         prefix="$2"
         shift 2
         ;;
      --disable-shared)
         shared_flag="--disable-shared"
         shift
         ;;
      *)
         echo "Usage: $0 [--prefix DIR] [--disable-shared]"
         exit 1
         ;;
   esac
done

script_dir=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$prefix"

# SPLASH uses the C HDF5 API only; static builds avoid runtime HDF5 deps.
# -fPIC is required when the final binary is linked with -pie (default PIELDFLAGS).
if [ -z "${CFLAGS:-}" ]; then
   export CFLAGS="-fPIC"
fi

"${script_dir}/install-pkg.sh" "$hdf5_url" "$prefix" "--enable-fortran=no $shared_flag"
