"""Utilities to iterate local files."""
import os
import sys
import tempfile
import urllib.parse
import urllib.request
from typing import Generator, Optional, Union

import pandas as pd

from access import decompress, download


def txt(url: str,
        file_from_zip_archive: str = "",
        buffering: int = 1024) -> Generator[str, None, None]:
    """
    Downloads, iterates and subsequently removes the file at a given URL.

    Args:
        url: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        buffering: The buffer size to process download, possible decompression
            and iteration of the file at.


    Yields:
        Lines of the file at a given URL.
    """
    if not os.path.exists(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            )):
        os.mkdir(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            ))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}",
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not os.path.exists(local_file_name):
        download.download_file(url, local_file_name, buffering)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress.decompress_gzip_file(
            local_file_name, buffering)
    elif file_name_extension == ".zip":
        local_file_name = decompress.decompress_zip_file(
            local_file_name, file_from_zip_archive)
    with open(local_file_name, buffering=buffering) as local_file:
        for line in local_file:
            yield line.rstrip("\n")

    os.remove(local_file_name)

    if not os.listdir(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            )):
        os.rmdir(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            ))


def tabular_txt(url: str,
                file_from_zip_archive: str = "",
                delimiter: str = "",
                header: int = 0,
                skiprows: int = 0,
                usecols: Optional[list[Union[int, str]]] = None,
                chunksize: int = 1024) -> Generator[pd.Series, None, None]:
    """
    Downloads, iterates and subsequently removes the tabular file at a URL.

    Args:
        url: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        delimiter: The delimiter separating values from different columns.
        header: The line number of the table header.
        skiprows: Number of initial rows of the file to skip before header.
        usecols: The table columns to consider.
        chunksize: The buffer size to process download, decompression and
            iteration of the file at.

    Yields:
        Rows of the file at a given URL.
    """
    if not os.path.exists(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            )):
        os.mkdir(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            ))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}",
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not os.path.exists(local_file_name):
        download.download_file(url, local_file_name, chunksize)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress.decompress_gzip_file(
            local_file_name, chunksize)
    elif file_name_extension == ".zip":
        local_file_name = decompress.decompress_zip_file(
            local_file_name, file_from_zip_archive)

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
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            )):
        os.rmdir(
            os.path.join(
                tempfile.gettempdir(),
                f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-{os.getpid()}"
            ))
