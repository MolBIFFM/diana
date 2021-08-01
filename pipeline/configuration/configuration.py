import os
import tempfile

CHUNK_SIZE = 2**16
USER_AGENT = "pipeline"
DOWNLOAD_DIRECTORY = os.path.join(tempfile.gettempdir(), "pipeline")
