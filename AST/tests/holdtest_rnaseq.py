import os
from pathlib import Path

import pytest

from AST.genelab_api import download_isa
import AST

def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(
        argnames, [[funcargs[name] for name in argnames] for funcargs in funcarglist]
    )

@pytest.fixture(scope='module')
def local_isa():
    yield {
        '194':Path(__file__).parent / Path('assets/GLDS-194_metadata_GLDS-194-ISA.zip') # RNASeq
        }

class TestClass:
    # a map specifying multiple argument sets for a test method
    params = {
        "test_from_repo_fetch_isazip": [
            dict(
              accession='GLDS-104',
              alternate_url=False,
              expected="GLDS-104_metadata_GLDS-104-ISA.zip"
              ),
            dict(
              accession='GLDS-194',
              alternate_url=True,
              expected="GLDS-194_metadata_GLDS-194-ISA.zip"
              ),
            dict(
              accession='GLDS-373',
              alternate_url=True,
              expected="GLDS-373_metadata_GLDS-373-ISA.zip"
              ),
        ],
        "test_from_local_isazip_load_isa": [
            dict(
              accession='194'
              ),
        ],
        "test_from_local_isazip_generate_proto_runsheet": [
            dict(
              accession='194',
              template=None,
              ),
        ],
    }

    def test_from_repo_fetch_isazip(self, accession, alternate_url, expected, tmpdir):
        print(os.chdir(tmpdir))
        isazip = download_isa(accession=accession, alternate_url = alternate_url)
        assert isazip == expected

    def test_from_local_isazip_load_isa(self, accession, local_isa):
        print(os.getcwd())
        proto_runsheet = AST.load_isazip( isazip = local_isa[accession] )
        #assert proto_runsheet == "GLDS-5_metadata_E-GEOD-4209-ISA.csv"

    def test_from_local_isazip_generate_proto_runsheet(self, accession, template, local_isa):
        print(os.getcwd())
        proto_runsheet = AST.isazip_to_runsheet( isazip = local_isa[accession], template = template)
        #assert proto_runsheet == "GLDS-5_metadata_E-GEOD-4209-ISA.csv"
