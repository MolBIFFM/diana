import gzip
import os
import tempfile
import urllib.parse
import urllib.request
import zipfile

import pandas as pd

def download_file(url):
    remote_file_name = os.path.split(urllib.parse.urlparse(url).path)[1]
    local_file_name = os.path.join(tempfile.gettempdir(), remote_file_name)

    if not os.path.exists(local_file_name):
        request = urllib.request.Request(
            url, headers={"User-Agent": "https://github.com/lucasfein/pipeline"})

        with urllib.request.urlopen(request) as response:
            with open(local_file_name, "wb") as local_file:
                chunk = response.read(1048576)
                while chunk:
                    local_file.write(chunk)
                    chunk = response.read(1048576)

    return local_file_name


def download_gzip_file(url):
    compressed_file_name = download_file(url)
    decompressed_file_name = os.path.splitext(compressed_file_name)[0]

    with gzip.open(compressed_file_name, "rb") as compressed_file:
        with open(decompressed_file_name, "wb") as decompressed_file:
            chunk = compressed_file.read(1048576)
            while chunk:
                decompressed_file.write(chunk)
                chunk = compressed_file.read(1048576)

    os.remove(compressed_file_name)

    return decompressed_file_name


def download_zip_file(url):
    compressed_file_name = download_file(url)

    with zipfile.ZipFile(compressed_file_name) as archive:
        decompressed_file_name = archive.extract(archive.namelist()[0],
                                                 path=tempfile.gettempdir())

    os.remove(compressed_file_name)

    return decompressed_file_name


def iterate_tabular_data(url, delimiter=None, header=None, usecols=[]):
    file_name_extension = os.path.splitext(urllib.parse.urlparse(url).path)[1]

    if file_name_extension == ".gz":
        local_file_name = download_gzip_file(url)
    elif file_name_extension == ".zip":
        local_file_name = download_zip_file(url)
    else:
        local_file_name = download_file(url)

    if os.path.splitext(local_file_name)[1] == ".csv":
        delimiter = ","
    elif os.path.splitext(local_file_name)[1] == ".tsv":
        delimiter = "\t"

    return pd.read_csv(local_file_name,
                       sep=delimiter,
                       header=header,
                       usecols=usecols).iterrows()
