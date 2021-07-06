from urllib.request import urlopen, quote, urlretrieve
from json import loads
import tempfile
import zipfile
from pathlib import Path

from isatools import isatab

def read_json(url):
    with urlopen(url) as response:
        return loads(response.read().decode())

def get_unzipped_isa_dir(isa_zip_path: str) -> str:
    """ Unzips ISA and places into a tmp contents folder.
    Returns path to temporary directory holding ISA zip file contents.

    Occasionally, the isa zip folder contains a nested directory of the same name.
    In these case, returns this nested folder containing isa files.

    :param isa_zip_path: path to isa zip file
    """
    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(isa_zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    # check if isa_temp_dir contains the files or an additional nested directory
    contents = list(Path(temp_dir).iterdir())
    if len(contents) == 1: # format for nested isa zip files
        temp_dir = contents[0] # use nest dir that actually contains the files
    return temp_dir

def extract_i_file(isa_zip_path: str) -> Path:
    """ Extracts i_file from an isa zip file.  Raises an exception if one and only one i_file is not found
    """
    temp_dir = get_unzipped_isa_dir(isa_zip_path)
    i_file = list(Path(temp_dir).glob("i_*"))
    assert len(i_file) == 1
    return i_file[0]


def load_isazip(isazip: Path) -> dict:
    """ For a given isa_zip_path, loads represenations of the data as follows:

    i_* file: a dictionary of keys:DataFrame
    s_* file: a single DataFrame
    a_* file(s): a list of DataFrames

    :param isa_zip_path: isa zip file path.  As generated at GeneLab.
    """
    isa_temp_dir = get_unzipped_isa_dir(isazip)
    # find files using glob
    i_file_glob = list(Path(isa_temp_dir).glob("i_*"))
    s_file_glob = list(Path(isa_temp_dir).glob("s_*"))
    a_file_glob = list(Path(isa_temp_dir).glob("a_*"))

    # ensure only one of each file was found in the unzipped file
    assert 1 == len(i_file_glob) == len(s_file_glob), f"Error: there should be one and only 1 i_*, and s_* file"

    # pull the single file from the list of one
    i_file = i_file_glob[0]
    s_file = s_file_glob[0]

    # setup dict
    df_dict = dict()

    # parse into tables (study and assay files) or dictionary of tables (investigation file)
    df_dict["investigation"] = isatab.load_investigation(fp=i_file.open())
    df_dict["studies"] = isatab.load_table(fp=s_file.open())
    df_dict["assays"] = [isatab.load_table(fp=a_file.open()) for a_file in a_file_glob]

    return df_dict
