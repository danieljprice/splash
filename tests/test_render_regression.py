"""CLI render regression: compare splash PNG outputs to expected perceptual hashes."""
from pathlib import Path
import json

from PIL import Image
import imagehash
import pytest


def load_expectations(testdata_root: Path) -> dict:
    path = testdata_root / "expected_hashes.json"
    assert path.is_file(), f"missing {path}"
    with path.open() as fh:
        data = json.load(fh)
    assert data.get("hash") == "phash"
    return data


def phash_distance(img_a: Image.Image, img_b: Image.Image) -> int:
    return imagehash.phash(img_a) - imagehash.phash(img_b)


def phash_of(path: Path) -> imagehash.ImageHash:
    return imagehash.phash(Image.open(path))


class TestRenderRegression:

    def test_expectations_file(self, splash_testdata):
        data = load_expectations(splash_testdata)
        assert "images" in data and data["images"]

    @pytest.mark.parametrize(
        "image_name",
        [
            "log_rho_render.png",
            "log_rho_sink0_render.png",
            "log_rho_sink1_render.png",
            "log_rho_sink2_render.png",
            "log_rho_v_render.png",
            "log_u_render.png",
        ],
    )
    def test_render_phash(self, image_name, splash_testdata, work_dir):
        data = load_expectations(splash_testdata)
        spec = data["images"][image_name]
        expected = imagehash.hex_to_hash(spec["phash"])
        max_distance = int(spec["max_distance"])

        target_path = work_dir / image_name
        assert target_path.is_file(), f"missing render output: {target_path}"

        got = phash_of(target_path)
        dist = got - expected
        assert dist <= max_distance, (
            f"{image_name}: phash distance {dist} > {max_distance} "
            f"(got {got}, expected {expected})"
        )

        # Negative control: must not look like the deliberately bad image
        neg = data["negative_control"]
        bad_path = splash_testdata / "control_images" / neg["file"]
        assert bad_path.is_file(), f"missing negative control {bad_path}"
        bad_hash = imagehash.hex_to_hash(neg["phash"])
        min_bad = int(neg["min_distance_vs_any_target"])
        bad_dist = got - bad_hash
        assert bad_dist > min_bad, (
            f"{image_name}: too similar to negative control "
            f"(distance {bad_dist} <= {min_bad})"
        )
