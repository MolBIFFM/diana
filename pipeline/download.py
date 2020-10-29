import gzip
import pathlib
import tempfile
import urllib.parse
import urllib.request
import zipfile

import pandas as pd

import pipeline.configuration


def download_file(url):
    file_name_suffix = "".join(
        pathlib.PurePosixPath(urllib.parse.urlparse(url).path).suffixes)

    request = urllib.request.Request(
        url, headers={"User-Agent": pipeline.configuration.USER_AGENT})
    with urllib.request.urlopen(request) as response:
        with tempfile.NamedTemporaryFile("wb",
                                         delete=False,
                                         suffix=file_name_suffix) as file:
            file.write(response.read())
    return pathlib.PurePosixPath(file.name)


def download_gzip_file(url):
    compressed_file_name = download_file(url)
    file_name_suffix = "".join(compressed_file_name.suffixes[:-1])

    with gzip.open(compressed_file_name, "rb") as compressed_file:
        with tempfile.NamedTemporaryFile("wb",
                                         delete=False,
                                         suffix=file_name_suffix) as file:
            file.write(compressed_file.read())
    return pathlib.PurePosixPath(file.name)


def download_zip_file(url):
    compressed_file_name = download_file(url)
    file_name_suffix = "".join(compressed_file_name.suffixes[:-1])

    with zipfile.ZipFile(compressed_file_name) as archive:
        with archive.open(archive.namelist()[0]) as compressed_file:
            with tempfile.NamedTemporaryFile("wb",
                                             delete=False,
                                             suffix=file_name_suffix) as file:
                file.write(compressed_file.read())
    return pathlib.PurePosixPath(file.name)


def iterate_tabular_data(url, delimiter=None, header=None, usecols=[]):
    file_name = pathlib.PurePosixPath(urllib.parse.urlparse(url).path)

    if file_name.suffix == ".gz":
        local_file_name = download_gzip_file(url)
    elif file_name.suffix == ".zip":
        local_file_name = download_zip_file(url)
    else:
        local_file_name = download_file(url)

    if local_file_name.suffix == ".csv":
        delimiter = ","
    elif local_file_name.suffix in (".tsv", ".tab3"):
        delimiter = "\t"

    return pd.read_csv(local_file_name,
                       sep=delimiter,
                       header=header,
                       usecols=usecols).iterrows()
