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
            url,
            headers={"User-Agent": "https://github.com/lucasfein/pipeline"})

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

    if not os.path.exists(decompressed_file_name):
        with gzip.open(compressed_file_name, "rb") as compressed_file:
            with open(decompressed_file_name, "wb") as decompressed_file:
                chunk = compressed_file.read(1048576)
                while chunk:
                    decompressed_file.write(chunk)
                    chunk = compressed_file.read(1048576)

    return decompressed_file_name


def download_zip_file(url, file=None):
    compressed_file_name = download_file(url)

    with zipfile.ZipFile(compressed_file_name) as archive:
        if not file:
            file = archive.namelist()[0]

        if not os.path.exists(os.path.join(tempfile.gettempdir(), file)):
            decompressed_file_name = archive.extract(
                file, path=tempfile.gettempdir())
        else:
            decompressed_file_name = os.path.join(tempfile.gettempdir(), file)

    return decompressed_file_name


def read_tabular_data(url, file=None, delimiter=None, header=None, usecols=[]):
    file_name = urllib.parse.urlparse(url).path
    file_name_extension = os.path.splitext(file_name)[1]

    if file_name_extension == ".gz":
        local_file_name = download_gzip_file(url)
    elif file_name_extension == ".zip":
        local_file_name = download_zip_file(url, file)
    else:
        local_file_name = download_file(url)

    if os.path.splitext(local_file_name)[1] == ".csv":
        delimiter = ","
    elif os.path.splitext(local_file_name)[1] == ".tsv":
        delimiter = "\t"

    return pd.read_csv(local_file_name,
                       sep=delimiter,
                       header=header,
                       usecols=usecols,
                       engine="c").iterrows()
