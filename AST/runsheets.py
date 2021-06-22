from pathlib import Path

from AST.utils import load_isazip

class RunSheet():
    def __init__(self):
        self.x = 1

class MicroarrayRunSheet(RunSheet):
    def __init__(self):
        self.x = 1

def isazip_to_runsheet(isazip: Path, template: str = None) -> RunSheet:
    """ Reads an isazip file and returns runsheets corresponding to each valid assay

    :param isazip: The isazip file to parse.
    :param template: A peppy template to use.  If None, a proto_runsheet is returned instead
    """
    dfs = load_isazip(isazip)

    # using isatools parser load files
