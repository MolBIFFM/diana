"""Utilities to iterate local files."""
import gzip
import os
import re
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import zipfile
from typing import Collection, Iterator, Optional

import pandas as pd

from access import decompress, download


def txt(file: str,
        file_from_zip_archive: Optional[re.Pattern[str]] = None,
        buffering: int = 1048576,
        pause: float = 60.0) -> Iterator[str]:
    """
    Downloads, iterates and subsequently removes the file at a given path.

    Args:
        file: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        buffering: The buffer size to process download, decompression
            and iteration of the file at.
        pause: The number of seconds to wait after a failed download.

    Yields:
        Lines of the file at a given path.
    """
    # Download the file from a URL to a temporary file location and decompress
    # the temporary file.
    if urllib.parse.urlparse(file).scheme in ("ftp", "http", "https"):
        if not os.path.exists(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}")):
            os.mkdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}"))

        local_file_name = os.path.join(
            tempfile.gettempdir(),
            f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
            f"{os.getpid()}",
            os.path.split(urllib.parse.urlparse(file).path)[1],
        )

        while True:
            try:
                if not os.path.exists(local_file_name):
                    download.download_file(file,
                                           local_file_name,
                                           buffering,
                                           pause=pause)

                file_name_extension = os.path.splitext(local_file_name)[1]

                if file_name_extension == ".gz":
                    local_file_name = decompress.decompress_gzip_file(
                        local_file_name, buffering)
                elif file_name_extension == ".zip":
                    local_file_name = decompress.decompress_zip_file(
                        local_file_name, file_from_zip_archive)
                break

            except (gzip.BadGzipFile, zipfile.BadZipFile):
                time.sleep(pause)

    # Decompress the local file.
    elif os.path.isfile(file):
        local_file_name = file

        file_name_extension = os.path.splitext(local_file_name)[1]

        if file_name_extension == ".gz":
            local_file_name = decompress.decompress_gzip_file(local_file_name,
                                                              buffering,
                                                              remove=False)
        elif file_name_extension == ".zip":
            local_file_name = decompress.decompress_zip_file(
                local_file_name, file_from_zip_archive, remove=False)

    else:
        return

    # Yield lines of the file.
    with open(local_file_name, buffering=buffering,
              encoding="utf-8") as local_file:
        for line in local_file:
            yield line.rstrip("\n")

    # Remove the downloaded file from the temporary file location.
    if urllib.parse.urlparse(file).scheme in ("ftp", "http", "https"):
        os.remove(local_file_name)

        if not os.listdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}")):
            os.rmdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}"))

    # Remove the file decompressed from a local compressed file.
    elif os.path.isfile(file) and not os.path.samefile(file, local_file_name):
        os.remove(local_file_name)


def tabular_txt(file: str,
                file_from_zip_archive: Optional[re.Pattern[str]] = None,
                delimiter: str = "",
                header: int = 0,
                skiprows: int = 0,
                usecols: Optional[Collection[int] | Collection[str]] = None,
                chunksize: int = 1048576,
                rows: int = 8192,
                pause: float = 60.0) -> Iterator[pd.Series]:
    """
    Downloads, iterates and subsequently removes the tabular file at a path.

    Args:
        file: The file location.
        file_from_zip_archive: The file from the zip archive to extract.
        delimiter: The delimiter separating values from different columns.
        header: The line number of the table header.
        skiprows: Number of initial rows of the file to skip before header.
        usecols: The table columns to consider.
        chunksize: The buffer size to process download and gzip decompression of
            the file with.
        rows: The number of rows to read at once from the table.
        pause: The number of seconds to wait after a failed download.

    Yields:
        Rows of the file at a given path.
    """
    # Download the file from a URL to a temporary file location and decompress
    # the temporary file.
    if urllib.parse.urlparse(file).scheme in ("ftp", "http", "https"):
        if not os.path.exists(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}"
                    f"-{os.getpid()}")):
            os.mkdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}"))

        local_file_name = os.path.join(
            tempfile.gettempdir(),
            f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
            f"{os.getpid()}",
            os.path.split(urllib.parse.urlparse(file).path)[1],
        )

        while True:
            try:
                if not os.path.exists(local_file_name):
                    download.download_file(file, local_file_name, chunksize)

                file_name_extension = os.path.splitext(local_file_name)[1]

                if file_name_extension == ".gz":
                    local_file_name = decompress.decompress_gzip_file(
                        local_file_name, chunksize)
                elif file_name_extension == ".zip":
                    local_file_name = decompress.decompress_zip_file(
                        local_file_name, file_from_zip_archive)
                break

            except (gzip.BadGzipFile, zipfile.BadZipFile):
                time.sleep(pause)

    # Decompress the local file.
    elif os.path.isfile(file):
        local_file_name = file

        file_name_extension = os.path.splitext(local_file_name)[1]

        if file_name_extension == ".gz":
            local_file_name = decompress.decompress_gzip_file(local_file_name,
                                                              chunksize,
                                                              remove=False)
        elif file_name_extension == ".zip":
            local_file_name = decompress.decompress_zip_file(
                local_file_name, file_from_zip_archive, remove=False)

    else:
        return

    # Determine the unspecified delimiter from the file extension.
    if not delimiter:
        if os.path.splitext(local_file_name)[1] == ".csv":
            delimiter = ","
        elif os.path.splitext(local_file_name)[1] == ".tsv":
            delimiter = "\t"

    # Yield rows of the file.
    for chunk in pd.read_csv(local_file_name,
                             delimiter=delimiter,
                             header=header,
                             skiprows=skiprows,
                             usecols=usecols,
                             chunksize=rows):
        for _, row in chunk.iterrows():
            yield row

    # Remove the downloaded file from the temporary file location.
    if urllib.parse.urlparse(file).scheme in ("ftp", "http", "https"):
        os.remove(local_file_name)

        if not os.listdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}")):
            os.rmdir(
                os.path.join(
                    tempfile.gettempdir(),
                    f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}-"
                    f"{os.getpid()}"))

    # Remove the file decompressed from a local compressed file.
    elif os.path.isfile(file) and not os.path.samefile(file, local_file_name):
        os.remove(local_file_name)
