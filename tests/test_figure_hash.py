from PIL import Image
import imagehash
from pathlib import Path
import os
import pytest


class TestImageHash():
    test_data_dir = Path(os.environ.get("SPLASH_DIR")) / 'data' 

    # Load your images
    img1 = Image.open(test_data_dir / 'out_2i.png')
    img2 = Image.open(test_data_dir / 'out_2i_2.png')

    # Generate perceptual hashes (pHash)
    hash1 = imagehash.phash(img1)
    hash2 = imagehash.phash(img2)

    print(f"Hash 1: {hash1}")
    print(f"Hash 2: {hash2}")

    # Calculate Hamming Distance (number of differing bits)
    # A distance of 0 means the images are perceptually identical.
    hamming_dist = hash1 - hash2
    print(f"Hamming Distance: {hamming_dist}")

    # Evaluate similarity based on a standard threshold
    threshold = 10
    assert hamming_dist <= threshold