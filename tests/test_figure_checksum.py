import hashlib

def get_file_sha256(filename):
    sha256_hash = hashlib.sha256()

    # Read the file in binary mode
    with open(filename, "rb") as f:
        # Read in 64KB chunks to efficiently handle large files
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)

    return sha256_hash.hexdigest()

get_file_sha256('out_2i.png')
