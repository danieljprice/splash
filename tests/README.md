# CLI render regression tests
# ===========================
#
# This directory contains both Fortran unit tests
# and a regression testing script that runs `splash` on a binary dump from phantom, 
# writes pngs, and compare `imagehash.phash` distances to expectations in
# https://github.com/danieljprice/splash-testdata (latest `main`).
#
# The Phantom dump used for renders lives in that repo (`datafiles/binary_00000`)
#
# ## Dependencies (CI uses apt + unpinned pip)
#
# Ubuntu / Debian:
#   sudo apt-get install -y python3-pytest python3-pil python3-pip
#   pip3 install imagehash
#
# macOS:
#   brew install python
#   pip3 install pytest pillow imagehash
#
# ## Workflow run
#
# Regression testing is performed automatically in github actions
# via .github/workflows/regresssion.yml
#  
# ## Local run
#
# 1. Build splash and put `bin/` on PATH (and giza libs on the library path).
# 2. Clone fixtures:
#
#        git clone https://github.com/danieljprice/splash-testdata.git
#        mkdir -p test_data
#
# 3. Generate renders and run pytest:
#
#        export SPLASH_TESTDATA=$PWD/splash-testdata
#        export SPLASH_WORK_DIR=$PWD/test_data
#        ./scripts/run_render_regression.sh
#        python3 -m pytest -v tests/
#
# ## Updating expected hashes
#
# Recompute on Linux is preferred. Update
# `control_images/` and `expected_hashes.json` in splash-testdata and push
# to that repository's `main`.
