import os
from pathlib import Path

import pytest

from AST import isa_cleanup
from AST.utils import get_unzipped_isa_dir, extract_i_file

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
        '5':Path(__file__).parent / Path('assets/GLDS-5_metadata_E-GEOD-4209-ISA.zip'),
        '121':Path(__file__).parent / Path('assets/BRIC-16_CYTO-ISA.zip'),
        }

class TestClass:
    LOCAL_ISA_ZIP_SET = [dict(accession='5'), dict(accession='121')]
    # a map specifying multiple argument sets for a test method
    params = {
        "test_extract_i_file": [
            dict(
                accession='5',
                expected_i_file='i_E-GEOD-4209_investigation.txt'
            ),
        ],
        "test_isa_cleanup_read_i_file": LOCAL_ISA_ZIP_SET,
        "test_isa_cleanup_remove_double_quotes": [dict(accession='121')],

    }
    def test_extract_i_file(self, accession, expected_i_file, local_isa):
        i_file = extract_i_file(local_isa[accession])
        assert i_file.name == expected_i_file

    def test_isa_cleanup_remove_double_quotes(self, accession, local_isa):
        i_file = extract_i_file(local_isa[accession])
        new_i_file = isa_cleanup.cleanup_double_quotes(i_file)
        assert new_i_file == "null"
