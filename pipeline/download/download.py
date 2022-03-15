import gzip
import os
import re
import tempfile
from typing import Generator, Union
import urllib.parse
import urllib.request
import sys
import zipfile

import pandas as pd


def download_file(url: str, local_file_name: str, size: int = 1048576) -> None:
    """
    Downloads a file from a URL.

    Args:
        url: The file location.
        local_file_name: The local file name to download the file to.
        size: The buffer size to process the download at.
    """
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": os.path.splitext(os.path.basename(sys.argv[0]))[0]
        })

    with urllib.request.urlopen(request) as response:
        with open(local_file_name, "wb") as local_file:
            while chunk := response.read(size):
                local_file.write(chunk)


def decompress_gzip_file(compressed_file_name: str, size: int = 1048576) -> str:
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
                    tempfile.gettempdir(), "{}-{}".format(
                        os.path.splitext(os.path.basename(sys.argv[0]))[0],
                        os.getpid()), file)):
            decompressed_file_name = archive.extract(
                file,
                path=os.path.join(
                    tempfile.gettempdir(), "{}-{}".format(
                        os.path.splitext(os.path.basename(sys.argv[0]))[0],
                        os.getpid())))
        else:
            decompressed_file_name = os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid()), file)

    os.remove(compressed_file_name)
    return decompressed_file_name


def txt(url: str,
        file_from_zip_archive: str = "",
        buffering: int = 1048576) -> Generator[str, None, None]:
    """
    Downloads, iterates and subsequently removes the file at a given URL.

    Args:
        url: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        buffering: The buffer size to process download, possible decompression and iteration of the file at.


    Yields:
        Lines of the file at a given URL.
    """
    if not os.path.exists(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid()))):
        os.mkdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid())))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        "{}-{}".format(
            os.path.splitext(os.path.basename(sys.argv[0]))[0], os.getpid()),
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not (os.path.exists(local_file_name)):
        download_file(url, local_file_name, buffering)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress_gzip_file(local_file_name, buffering)
    elif file_name_extension == ".zip":
        local_file_name = decompress_zip_file(local_file_name,
                                              file_from_zip_archive)
    with open(local_file_name, buffering=buffering) as local_file:
        for line in local_file:
            yield line.rstrip("\n")

    os.remove(local_file_name)

    if not os.listdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid()))):
        os.rmdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid())))


def tabular_txt(url: str,
                file_from_zip_archive: str = "",
                delimiter: str = "",
                header: int = 0,
                skiprows: int = 0,
                usecols: list[Union[int, str]] = [],
                chunksize: int = 1048576) -> Generator[pd.Series, None, None]:
    """
    Downloads, iterates and subsequently removes the tabular file at a given URL.

    Args:
        url: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        delimiter: The delimiter separating values from different columns.
        header: The line number of the table header.
        skiprows: Number of initial rows of the file to skip before header.
        usecols: The table columns to consider.
        chunksize: The buffer size to process download, decompression and iteration of the file at.


    Yields:
        Rows of the file at a given URL.
    """
    if not os.path.exists(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid()))):
        os.mkdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid())))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        "{}-{}".format(
            os.path.splitext(os.path.basename(sys.argv[0]))[0], os.getpid()),
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not (os.path.exists(local_file_name)):
        download_file(url, local_file_name, chunksize)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress_gzip_file(local_file_name, chunksize)
    elif file_name_extension == ".zip":
        local_file_name = decompress_zip_file(local_file_name,
                                              file_from_zip_archive)

    if not delimiter:
        if os.path.splitext(local_file_name)[1] == ".csv":
            delimiter = ","
        elif os.path.splitext(local_file_name)[1] == ".tsv":
            delimiter = "\t"

    for chunk in pd.read_csv(
            local_file_name,
            delimiter=delimiter,
            header=header,
            skiprows=skiprows,
            usecols=usecols,
            chunksize=chunksize,
    ):
        for _, row in chunk.iterrows():
            yield row

    os.remove(local_file_name)

    if not os.listdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid()))):
        os.rmdir(
            os.path.join(
                tempfile.gettempdir(), "{}-{}".format(
                    os.path.splitext(os.path.basename(sys.argv[0]))[0],
                    os.getpid())))
