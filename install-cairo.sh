#!/bin/bash
#
# Script for splash 2.x that retrieves and installs
# both cairo and pixman
#
# Usage: install-cairo.sh [--prefix DIR] [--disable-shared]
#
set -e

pixman_url="https://cairographics.org/releases/pixman-0.46.4.tar.gz"
cairo_url="https://cairographics.org/releases/cairo-1.18.4.tar.xz"
prefix="$PWD/giza"
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

mkdir -p "$prefix"

# pixman >= 0.43 and cairo >= 1.17 use meson; 1.18.4 requires meson >= 1.3.0
# meson >= 1.3 requires Python >= 3.7; manylinux system python3 is often 3.6
meson_version_ok() {
   command -v meson >/dev/null 2>&1 || return 1
   meson --version 2>/dev/null | awk -F. \
      'NR==1 { exit !($1>1 || ($1==1 && $2>=3)) } END { if (NR==0) exit 1 }'
}

meson_python() {
   if [ -n "$MESON_PYTHON" ] && [ -x "$MESON_PYTHON" ]; then
      echo "$MESON_PYTHON"
      return 0
   fi
   local py
   for py in /opt/python/cp311-cp311/bin/python /opt/python/cp310-cp310/bin/python \
             /opt/python/cp39-cp39/bin/python python3; do
      if [ -x "$py" ] && "$py" -c 'import sys; sys.exit(sys.version_info < (3, 7))' 2>/dev/null; then
         echo "$py"
         return 0
      fi
   done
   echo python3
}

ensure_meson() {
   local py pip_args pip_user_bin py_bin
   py=$(meson_python)
   py_bin="$(dirname "$py")"
   pip_user_bin="$("$py" -m site --user-base)/bin"
   export PATH="$py_bin:$pip_user_bin:$HOME/.local/bin:$PATH"
   if command -v brew >/dev/null 2>&1; then
      export PATH="$(brew --prefix)/bin:$PATH"
   fi
   if meson_version_ok; then
      return 0
   fi
   if [ "$(uname)" = Darwin ] && command -v brew >/dev/null 2>&1; then
      echo ":: installing meson via brew"
      brew install meson
      export PATH="$(brew --prefix)/bin:$PATH"
      if meson_version_ok; then
         return 0
      fi
   fi
   echo ":: installing meson >= 1.3.0 via pip ($py)"
   pip_args=()
   if [ "$(uname)" = Darwin ]; then
      # Homebrew Python is PEP 668 externally-managed
      pip_args+=(--break-system-packages)
   fi
   if [ "$(id -u)" -ne 0 ]; then
      pip_args+=(--user)
   fi
   "$py" -m pip install "${pip_args[@]}" 'meson>=1.3.0'
   export PATH="$py_bin:$pip_user_bin:$HOME/.local/bin:$PATH"
   if ! meson_version_ok; then
      echo "ERROR: meson >= 1.3.0 required but got: $(command -v meson 2>/dev/null || echo 'not found')"
      meson --version 2>&1 || true
      exit 1
   fi
}

download_file() {
   local url="$1"
   local distfile
   distfile=$(basename "$url")
   if [ ! -f "$distfile" ]; then
      if type -p wget > /dev/null 2>&1; then
         wget "$url"
      else
         curl -LO "$url"
      fi
   fi
   if [ ! -f "$distfile" ]; then
      echo "ERROR: failed to download $url"
      exit 1
   fi
   echo "$distfile"
}

ensure_meson
MESON_PY=$(meson_python)
export PATH="$(dirname "$MESON_PY"):$( "$MESON_PY" -m site --user-base)/bin:$HOME/.local/bin:$PATH"
MESON=meson

if [ -n "$shared_flag" ]; then
   libtype=static
else
   libtype=both
fi

# --- pixman (meson) ---
pixman_dist=$(download_file "$pixman_url")
pixman_dir="${pixman_dist/.tar.gz/}"
if [ ! -d "$pixman_dir" ]; then
   tar xzf "$pixman_dist"
fi
echo ":: installing pixman from $pixman_dir"
cd "$pixman_dir"
# Force lib/ (not lib64/) so consumers can use a stable path on RHEL/manylinux
"$MESON" setup build --prefix="$prefix" --libdir=lib --default-library="$libtype" -Dgtk=disabled
ninja -C build
ninja -C build install
cd ..

export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

# --- cairo (meson); giza uses Xlib, not xcb ---
cairo_dist=$(download_file "$cairo_url")
cairo_dir="${cairo_dist/.tar.xz/}"
if [ ! -d "$cairo_dir" ]; then
   tar -Jxf "$cairo_dist"
fi
echo ":: installing cairo from $cairo_dir"
cd "$cairo_dir"
"$MESON" setup build --prefix="$prefix" --libdir=lib --default-library="$libtype" \
   -Dxlib=enabled -Dxcb=disabled -Dgtk2-utils=disabled
ninja -C build
ninja -C build install
cd ..

# Resolve where static libs landed (lib vs lib64)
cairo_libdir=""
for d in "$prefix/lib" "$prefix/lib64"; do
   if [ -f "$d/libcairo.a" ] && [ -f "$d/libpixman-1.a" ]; then
      cairo_libdir="$d"
      break
   fi
done

test -f "$prefix/include/cairo/cairo.h" || { echo "ERROR: cairo header missing after install"; exit 1; }
if [ -z "$cairo_libdir" ]; then
   echo "ERROR: libcairo.a / libpixman-1.a missing under $prefix/lib or $prefix/lib64"
   exit 1
fi
echo "CAIRO_LIBDIR=$cairo_libdir"

echo "export PKG_CONFIG_PATH=\"$cairo_libdir/pkgconfig:\${PKG_CONFIG_PATH:-}\""
echo "type \"make\" to compile SPLASH"; echo;
