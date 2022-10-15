"""Utilities to download files."""
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


def download_file(url: str,
                  local_file_name: str,
                  size: int = 8192,
                  pause: float = 60.0) -> None:
    """
    Downloads a file from a URL, reattempting failed or timed out downloads.

    Args:
        url: The file location.
        local_file_name: The local file name to download the file to.
        size: The buffer size to process the download at.
        pause: The number of seconds to wait before reattempting a failed or
            timed out download.
    """
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": os.path.splitext(os.path.basename(sys.argv[0]))[0]
        })

    while True:
        try:
            with urllib.request.urlopen(
                    request, context=ssl.create_default_context()) as response:
                with open(local_file_name, "wb") as local_file:
                    while chunk := response.read(size):
                        local_file.write(chunk)
            break

        except (urllib.error.URLError, TimeoutError):
            time.sleep(pause)
