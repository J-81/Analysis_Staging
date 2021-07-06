from urllib.request import quote
from re import search
from pathlib import Path
import requests

from AST.utils import read_json

# Function to pull metadata zip from GeneLab
# Credit to Kirill Grigorev
GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
GLDS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/data/"
FILELISTINGS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/filelistings/"
ISA_ZIP_REGEX = r'.*_metadata_.*[_-]ISA\.zip$'

def get_isa(accession: str):
    """ Returns isa filename as well as GeneLab URLS from the associated file listing

    :param accession: GLDS accession ID, e.g. GLDS-194
    """
    glds_json = read_json(GLDS_URL_PREFIX + accession)
    try:
        _id = glds_json[0]["_id"]
        print(f"File Listing URL: {FILELISTINGS_URL_PREFIX + _id}")
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    isa_entries = [
        entry for entry in read_json(FILELISTINGS_URL_PREFIX + _id)
        if search(ISA_ZIP_REGEX, entry["file_name"])
    ]
    if len(isa_entries) == 0:
        raise ValueError("Unexpected: no ISAs found")
    elif len(isa_entries) > 1:
        raise ValueError("Unexpected: multiple files match the ISA regex")
    else:
        entry = isa_entries[0]
        version = entry["version"]
        url = GENELAB_ROOT + entry["remote_url"] + "?version={}".format(version)
        alt_url = (
            GENELAB_ROOT + "/genelab/static/media/dataset/" +
            quote(entry["file_name"]) + "?version={}".format(version)
        )
        return entry["file_name"], version, url, alt_url


def download_isa(accession: str, alternate_url: bool = False):
    """ Downloads isa for given accession number.

    :param accession: GLDS accession number, e.g. GLDS-194
    :param alternate_url: if true, uses alternative url, both alternate and default url should fetch the same file
    """
    print(f"Accessing GeneLab API for ISA file. Accesion: {accession}")
    filename ,_, url, alt_url  = get_isa(accession)
    if not Path(filename).is_file():
        print(f"Successfully retrieved ISA file location from API.")
        use_url = url if not alternate_url else alt_url
        if not alternate_url:
            print("WARNING: The default URL did not work in tests.  If it still fails use the alternate url!")
        print(f"Downloading from {use_url}.")
        r = requests.get(use_url)
        # If the response was successful, no Exception will be raised
        r.raise_for_status()
        with open(filename, "wb") as f:
            f.write(r.content)
        print(f"Finished downloading ISA file: {filename}")
    else:
        print(f"Already downloaded {filename} to current folder")
    return filename
