import gzip
import os
import re
import tempfile
import urllib.parse
import urllib.request
import uuid
import zipfile

import pandas as pd

SUBDIRECTORY_PREFIX = "pipeline-{}".format(uuid.uuid4())


def download_file(url, local_file_name):
    request = urllib.request.Request(url, headers={"User-Agent": "pipeline"})

    with urllib.request.urlopen(request) as response:
        with open(local_file_name, "wb") as local_file:
            while chunk := response.read(1024):
                local_file.write(chunk)


def decompress_gzip_file(compressed_file_name):
    decompressed_file_name = os.path.splitext(compressed_file_name)[0]

    if not os.path.exists(decompressed_file_name):
        with gzip.open(compressed_file_name, "rb") as compressed_file:
            with open(decompressed_file_name, "wb") as decompressed_file:
                while chunk := compressed_file.read(1024):
                    decompressed_file.write(chunk)

    os.remove(compressed_file_name)
    return decompressed_file_name


def decompress_zip_file(compressed_file_name, file_from_zip_archive=None):
    with zipfile.ZipFile(compressed_file_name) as archive:
        if not file_from_zip_archive:
            file = archive.namelist()[0]
        else:
            regex = re.compile(file_from_zip_archive)
            file = next(filter(regex.match, archive.namelist()))

        if not os.path.exists(
                os.path.join(tempfile.gettempdir(), "{}-{}".format(
                    SUBDIRECTORY_PREFIX, os.getpid()), file)):
            decompressed_file_name = archive.extract(
                file,
                path=os.path.join(
                    tempfile.gettempdir(),
                    "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid())))
        else:
            decompressed_file_name = os.path.join(
                tempfile.gettempdir(), "{}-{}".format(SUBDIRECTORY_PREFIX,
                                                      os.getpid()), file)

    os.remove(compressed_file_name)
    return decompressed_file_name


def txt(url, file_from_zip_archive=None):
    if not os.path.exists(
            os.path.join(tempfile.gettempdir(), "{}-{}".format(
                SUBDIRECTORY_PREFIX, os.getpid()))):
        os.mkdir(
            os.path.join(tempfile.gettempdir(),
                         "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid())))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid()),
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not (os.path.exists(local_file_name)):
        download_file(url, local_file_name)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress_gzip_file(local_file_name)
    elif file_name_extension == ".zip":
        local_file_name = decompress_zip_file(local_file_name,
                                              file_from_zip_archive)

    with open(local_file_name, buffering=1024) as local_file:
        for line in local_file:
            yield line.rstrip("\n")

    os.remove(local_file_name)

    if not os.listdir(
            os.path.join(tempfile.gettempdir(), "{}-{}".format(
                SUBDIRECTORY_PREFIX, os.getpid()))):
        os.rmdir(
            os.path.join(tempfile.gettempdir(),
                         "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid())))


def tabular_txt(url,
                file_from_zip_archive=None,
                delimiter=None,
                header=None,
                skiprows=0,
                usecols=[]):
    if not os.path.exists(
            os.path.join(tempfile.gettempdir(), "{}-{}".format(
                SUBDIRECTORY_PREFIX, os.getpid()))):
        os.mkdir(
            os.path.join(tempfile.gettempdir(),
                         "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid())))

    local_file_name = os.path.join(
        tempfile.gettempdir(),
        "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid()),
        os.path.split(urllib.parse.urlparse(url).path)[1],
    )

    if not (os.path.exists(local_file_name)):
        download_file(url, local_file_name)

    file_name_extension = os.path.splitext(local_file_name)[1]
    if file_name_extension == ".gz":
        local_file_name = decompress_gzip_file(local_file_name)
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
            chunksize=1024,
    ):
        for _, row in chunk.iterrows():
            yield row

    os.remove(local_file_name)

    if not os.listdir(
            os.path.join(tempfile.gettempdir(), "{}-{}".format(
                SUBDIRECTORY_PREFIX, os.getpid()))):
        os.rmdir(
            os.path.join(tempfile.gettempdir(),
                         "{}-{}".format(SUBDIRECTORY_PREFIX, os.getpid())))
