"""Utilities to decompress local files."""
import gzip
import os
import re
import sys
import tempfile
import zipfile
from typing import Optional


def decompress_gzip_file(compressed_file_name: str, size: int = 8192) -> str:
    """
    Decompresses a gzip compressed file and removes the compressed file.

    Args:
        compressed_file_name: The file name of the compressed file.
        size: The buffer size to process the decompression at.

    Returns:
        The file name of the decompressed file.
    """
    decompressed_file_name = os.path.splitext(compressed_file_name)[0]

    if not os.path.exists(decompressed_file_name):
        with gzip.open(compressed_file_name, "rb") as compressed_file:
            with open(decompressed_file_name,
                      "wb",
                      buffering=size,
                      encoding="utf-8") as decompressed_file:
                while chunk := compressed_file.read(size):
                    decompressed_file.write(chunk)

    os.remove(compressed_file_name)
    return decompressed_file_name


def decompress_zip_file(
        compressed_file_name: str,
        file_from_zip_archive: Optional[re.Pattern] = None) -> str:
    """
    Decompresses a zip compressed file and removes the compressed file.

    Args:
        compressed_file_name: The file name of the compressed file.
        file_from_zip_archive: The file from the zip archive to extract. If
            this file is not specified or found, the first in the archive is
            used. The first file matching the regular expression is used.

    Returns:
        The file name of the decompressed file.
    """
    with zipfile.ZipFile(compressed_file_name) as archive:
        if not file_from_zip_archive:
            file = archive.namelist()[0]
        else:
            for name in archive.namelist():
                if file_from_zip_archive.fullmatch(name):
                    file = name
                    break
            else:
                file = archive.namelist()[0]

        if not os.path.exists(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}", file)):
            decompressed_file_name = archive.extract(
                file,
                path=os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}"))
        else:
            decompressed_file_name = os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                f"{os.getpid()}", file)

    os.remove(compressed_file_name)
    return decompressed_file_name
