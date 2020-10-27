import gzip
import pathlib
import tempfile
import urllib.parse
import urllib.request

import pandas as pd

import ppi.configuration


def download_file(url):
    request = urllib.request.Request(
        url, headers={"User-Agent": ppi.configuration.USER_AGENT})
    with urllib.request.urlopen(request) as response:
        with tempfile.NamedTemporaryFile("wb", delete=False) as file:
            file.write(response.read())
    return file.name


def download_gzip_file(url):
    with gzip.open(download_file(url), "rb") as content:
        with tempfile.NamedTemporaryFile("wb", delete=False) as file:
            file.write(content.read())
    return file.name


def get_tabular_data(url, header=None, usecols=[]):
    file_name = pathlib.PurePosixPath(urllib.parse.urlparse(url).path)

    if file_name.suffix == ".gz":
        local_file_name = download_gzip_file(url)
        compression = True
    else:
        local_file_name = download_file(url)
        compression = False

    if compression:
        delimiter_extension = file_name.suffixes[-2]
    else:
        delimiter_extension = file_name.suffixes[-1]

    if delimiter_extension == ".csv":
        delimiter = ","
    elif delimiter_extension == ".tsv":
        delimiter = "\t"
    else:
        delimiter = " "

    return pd.read_csv(local_file_name,
                       sep=delimiter,
                       engine="c",
                       header=header,
                       usecols=usecols,
                       memory_map=True).iterrows()
