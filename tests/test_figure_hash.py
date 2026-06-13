from PIL import Image
import imagehash
from pathlib import Path
import os
import pytest
from enum import IntEnum
import numpy as np

class Threshold(IntEnum):
    LOW    = 10
    MEDIUM =  7
    HIGH   =  4


SPLASH_DIR  = Path(os.environ.get("SPLASH_DIR"))
TARGET_DIR  = Path("./")
CONTROL_DIR = SPLASH_DIR / "data/control_images"
BAD_PLOT    = CONTROL_DIR / "bad_plot.png"



def hamming(img1, img2):
    h1 = imagehash.phash(img1).hash
    h2 = imagehash.phash(img2).hash
    return int(np.count_nonzero(h1 != h2))


@pytest.fixture
def image_pair():
    def _load(image_name):
        control = Image.open(CONTROL_DIR / image_name)
        target  = Image.open(TARGET_DIR  / image_name)
        return control, target
    return _load


class TestImageHash:

    def test_log_rho_render(self, image_pair):
        assert (TARGET_DIR / "log_rho_render.png").exists()
        control, target = image_pair("log_rho_render.png")
        assert hamming(control, target) <= Threshold.LOW
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW

    def test_log_rho_sink0_render(self, image_pair):
        assert (TARGET_DIR / "log_rho_sink0_render.png").exists()
        control, target = image_pair("log_rho_sink0_render.png")
        assert hamming(control, target) <= Threshold.LOW
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW

    def test_log_rho_sink1_render(self, image_pair):
        assert (TARGET_DIR / "log_rho_sink1_render.png").exists()
        control, target = image_pair("log_rho_sink1_render.png")
        assert hamming(control, target) <= Threshold.MEDIUM
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW

    def test_log_rho_sink2_render(self, image_pair):
        assert (TARGET_DIR / "log_rho_sink2_render.png").exists()
        control, target = image_pair("log_rho_sink2_render.png")
        assert hamming(control, target) <= Threshold.MEDIUM
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW

    def test_log_rho_v_render(self, image_pair):
        assert (TARGET_DIR / "log_rho_v_render.png").exists()
        control, target = image_pair("log_rho_v_render.png")
        assert hamming(control, target) <= Threshold.HIGH
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW

    def test_log_u_render(self, image_pair):
        assert (TARGET_DIR / "log_u_render.png").exists()
        control, target = image_pair("log_u_render.png")
        assert hamming(control, target) <= Threshold.LOW
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW