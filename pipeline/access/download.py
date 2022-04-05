"""Utilities to download files."""
import os
import sys
import time
import urllib.parse
import urllib.request


def download_file(url: str, local_file_name: str, size: int = 1024) -> None:
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

    while True:
        try:
            with urllib.request.urlopen(request) as response:
                with open(local_file_name, "wb") as local_file:
                    while chunk := response.read(size):
                        local_file.write(chunk)
            break

        except urllib.error.URLError:
            time.sleep(60)
