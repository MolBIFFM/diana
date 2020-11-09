import gzip
import os
import tempfile
import urllib.parse
import urllib.request
import zipfile

import pandas as pd


def download_file(url, local_file_name):
    request = urllib.request.Request(
        url, headers={"User-Agent": "lucas.fein@stud.uni-frankfurt.de"})

    with urllib.request.urlopen(request) as response:
        with open(local_file_name, "wb") as local_file:
            while chunk := response.read(1048576):
                local_file.write(chunk)


def decompress_gzip_file(compressed_file_name):
    decompressed_file_name = os.path.splitext(compressed_file_name)[0]

    if not os.path.exists(decompressed_file_name):
        with gzip.open(compressed_file_name, "rb") as compressed_file:
            with open(decompressed_file_name, "wb") as decompressed_file:
                while chunk := compressed_file.read(1048576):
                    decompressed_file.write(chunk)

    os.remove(compressed_file_name)
    return decompressed_file_name


def decompress_zip_file(compressed_file_name, file=None):
    with zipfile.ZipFile(compressed_file_name) as archive:
        if not file:
            file = archive.namelist()[0]

        if not os.path.exists(os.path.join(tempfile.gettempdir(), file)):
            decompressed_file_name = archive.extract(
                file, path=tempfile.gettempdir())
        else:
            decompressed_file_name = os.path.join(tempfile.gettempdir(), file)

    os.remove(compressed_file_name)
    return decompressed_file_name


def read_tabular_data(url, file=None, delimiter=None, header=None, usecols=[]): 
    if file:
        file_name_extension = os.path.splitext(file)[1]
        local_file_name = os.path.join(tempfile.gettempdir(), file)
    else:
        file_name = os.path.split(urllib.parse.urlparse(url).path)[1]
        file_name_extension = os.path.splitext(file_name)[1]
        local_file_name = os.path.join(tempfile.gettempdir(), file_name)

    if not os.path.exists(local_file_name):
        download_file(url, local_file_name)

    if file_name_extension == ".gz":
        local_file_name = decompress_gzip_file(local_file_name)
    elif file_name_extension == ".zip":
        local_file_name = decompress_zip_file(local_file_name, file)

    if os.path.splitext(local_file_name)[1] == ".csv":
        delimiter = ","
    elif os.path.splitext(local_file_name)[1] == ".tsv":
        delimiter = "\t"

    for chunk in pd.read_csv(local_file_name,
                       sep=delimiter,
                       header=header,
                       usecols=usecols,
                       chunksize=1048576):
        for _, row in chunk.iterrows():
            yield row
