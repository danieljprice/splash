from PIL import Image
import imagehash
from pathlib import Path
import os
import pytest
from enum import IntEnum
import numpy as np

class Threshold(IntEnum):
    HIGH     =  4  # Nearly identical, tight tolerance
    MEDIUM   =  7  # Minor platform differences expected
    LOW      = 10  # Moderate differences acceptable
    RELAXED  = 15  # Known sensitivity to platform/version differences


SPLASH_DIR  = Path(os.environ.get("SPLASH_DIR", Path(__file__).resolve().parents[1]))
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

    @pytest.mark.parametrize(
        "image_name, threshold",
        [
            ("log_rho_render.png", Threshold.LOW),
            ("log_rho_sink0_render.png", Threshold.LOW),
            ("log_rho_sink1_render.png", Threshold.RELAXED),
            ("log_rho_sink2_render.png", Threshold.RELAXED),
            ("log_rho_v_render.png", Threshold.LOW),
            ("log_u_render.png", Threshold.LOW),
        ]
    )
    def test_renders(self, image_name, threshold, image_pair):
        target_path = TARGET_DIR / image_name
        assert target_path.exists()
        control, target = image_pair(image_name)
        assert hamming(control, target) <= threshold
        assert hamming(Image.open(BAD_PLOT), target) > Threshold.LOW