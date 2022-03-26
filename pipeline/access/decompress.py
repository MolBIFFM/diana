"""Utilities to decompress local files."""
import gzip
import os
import re
import sys
import tempfile
import zipfile


def decompress_gzip_file(compressed_file_name: str, size: int = 4096) -> str:
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
            with open(decompressed_file_name, "wb") as decompressed_file:
                while chunk := compressed_file.read(size):
                    decompressed_file.write(chunk)

    os.remove(compressed_file_name)
    return decompressed_file_name


def decompress_zip_file(compressed_file_name: str,
                        file_from_zip_archive: str = "") -> str:
    """
    Decompresses a zip compressed file and removes the compressed file.

    Args:
        compressed_file_name: The file name of the compressed file.
        file_from_zip_archive: The file from the zip archive to extract.

    Returns:
        The file name of the decompressed file.
    """
    with zipfile.ZipFile(compressed_file_name) as archive:
        if not file_from_zip_archive:
            file = archive.namelist()[0]
        else:
            regex = re.compile(file_from_zip_archive)
            file = next(filter(regex.match, archive.namelist()))

        if not os.path.exists(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}",
                    file)):
            decompressed_file_name = archive.extract(
                file,
                path=os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
                ))
        else:
            decompressed_file_name = os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}",
                file)

    os.remove(compressed_file_name)
    return decompressed_file_name
