import os
import tempfile

CHUNK_SIZE = 1 * 2 ** 20
USER_AGENT = "pipeline"
DOWNLOAD_DIRECTORY = os.path.join(tempfile.gettempdir(), "pipeline")
