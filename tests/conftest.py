# Shared fixtures for splash black-box / CLI regression tests.
from pathlib import Path
import os

import pytest


@pytest.fixture(scope="session")
def splash_testdata():
    """Root of a splash-testdata checkout (env SPLASH_TESTDATA)."""
    env = os.environ.get("SPLASH_TESTDATA")
    if not env:
        pytest.skip("SPLASH_TESTDATA is not set")
    root = Path(env)
    if not root.is_dir():
        pytest.fail(f"SPLASH_TESTDATA does not exist: {root}")
    return root


@pytest.fixture(scope="session")
def work_dir():
    """Directory containing PNG outputs from run_render_regression.sh."""
    env = os.environ.get("SPLASH_WORK_DIR", ".")
    return Path(env)
