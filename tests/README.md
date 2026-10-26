# CLI regression tests

This directory contains Fortran unit tests plus pytest CLI tests that run
`splash` against dumps in https://github.com/danieljprice/splash-testdata
(latest `main`):

- `test_render_regression.py` compares PNG `imagehash.phash` distances
- `test_data_formats.py` runs `splash --labels` / `--header` (and temp
  `splash to ascii`) for one small dump per reader

The Phantom dump used for renders lives in that repo (`dataformats/phantom/binary_00000`).
Format fixtures live under `dataformats/` with expectations in
`expected_formats.json`.

## Dependencies (CI uses apt + unpinned pip)

Ubuntu / Debian:

```shell
sudo apt-get install -y python3-pytest python3-pil python3-pip
pip3 install imagehash
```

macOS:

```shell
brew install python
pip3 install pytest pillow imagehash
```

## Workflow run

Tests are run automatically in GitHub Actions via
`.github/workflows/test.yml`.

## Local run

1. Build splash and put `bin/` on PATH (and giza libs on the library path).
2. Clone fixtures:

```shell
git clone https://github.com/danieljprice/splash-testdata.git
mkdir -p test_data
```

3. Generate renders and run pytest:

```shell
export PATH=$PWD/bin:$PATH
export SPLASH_TESTDATA=$PWD/splash-testdata
export SPLASH_WORK_DIR=$PWD/test_data
./scripts/run_render_regression.sh
python3 -m pytest -v tests/
```

Format tests (`test_data_formats.py`) only need `splash` on PATH and
`SPLASH_TESTDATA`. They do not require the render script.

## Updating expected hashes

Update `control_images/` and `expected_hashes.json` in splash-testdata and push to that repository's `main`.
