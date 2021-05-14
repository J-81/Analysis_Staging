#! /usr/bin/env python
"""
WARNING: Testing with GLDS-194 implies the default url no longer works.
Alternative url does work!
"""
from __future__ import annotations
from urllib.request import urlopen, quote, urlretrieve
from json import loads
from re import search
import re
import argparse
import zipfile
import tempfile
import os
from pathlib import Path
from collections import defaultdict
import sys
import shutil
import importlib.resources
import traceback

import pandas as pd
from pandas import DataFrame
from isatools import isatab
from isatools.io import isatab_parser
from isatools.io.isatab_parser import ISATabRecord, NodeRecord
import requests
import peppy

################################################################################
# UTIILTIES
################################################################################

def _parse_args():
  """ Parse command line args.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('--accession', metavar='GLDS-001', required=True,
                      help='GLDS accesion number')
  parser.add_argument('--alternate-url', action="store_true", default=False,
                      help='Use alternate url, fetched by api script')
  parser.add_argument('--allow-missing-columns', action="store_true", default=False,
                      help='Creates protorun sheet even if missing metadata columns')
  #parser.add_argument('--to-runsheet', action="store_true", default=False,
  #                      help='Creates Runsheet based on ISA file and GeneLab API.  This mode autodetects the assay types')
  parser.add_argument('--to-RNASeq-runsheet', action="store_true", default=False,
                      help='Creates RNASeq Runsheet based on ISA file and GeneLab API')
  parser.add_argument('--alt-peppy-template', default=False,
                      help='Use alternative template')
  parser.add_argument('--to-Microarray-runsheet', action="store_true", default=False,
                      help='Creates Microarray Runsheet based on ISA file and GeneLab API')

  args = parser.parse_args()
  return args


def read_json(url):
    with urlopen(url) as response:
        return loads(response.read().decode())

################################################################################
# CONSTANTS
################################################################################

# Function to pull metadata zip from GeneLab
# Credit to Kirill Grigorev
GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
GLDS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/data/"
FILELISTINGS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/filelistings/"
ISA_ZIP_REGEX = r'.*_metadata_.*[_-]ISA\.zip$'

# assays with implemented runsheet templates
# tuple of measurement type and technology type
IMPLEMENTED_ASSAYS = [
    ("transcription profiling","RNA Sequencing (RNA-Seq)"),
    ("transcription profiling","Microarray"),
]

def get_assay_type(isa_zip_path: Path) -> (str, str):
    """ Returns the assay type for a given isa_zip_path

    :param isa_zip_path: Zip file containing genelab ISA files (handles nested directories as well as directly zipped files)
    """
    # unzip isa file
    temp_dir = _unzip_ISA(isa_zip_path)

    # find i_* file
    i_file_glob = list(Path(isa_temp_dir).glob("i_*"))
    assert 1 == len(i_file_glob), f"Error: there should be one and only 1 i_* file"
    i_file = i_file_glob[0]

    # get assay type
    assert len(i_df_dict["s_assays"]) == 1, "Should be length of one (dataframe)"
    i_assay_dict = query_assay_from_i_df_assays(
                                    i_df_dict["s_assays"][0],
                                    assay_measurement_type = "transcription profiling",
                                    assay_technology_type = "DNA Microarray"
                                    )

def load_isa(isa_zip_path: Path) -> (dict, DataFrame, list[DataFrame]):
    """ For a given isa_zip_path, loads represenations of the data as follows:

    i_* file: a dictionary of keys:DataFrame
    s_* file: a single DataFrame
    a_* file(s): a list of DataFrames

    :param isa_zip_path: isa zip file path.  As generated at GeneLab.
    """
    isa_temp_dir = _unzip_ISA(isazip)
    # find files using glob
    i_file_glob = list(Path(isa_temp_dir).glob("i_*"))
    s_file_glob = list(Path(isa_temp_dir).glob("s_*"))
    a_file_glob = list(Path(isa_temp_dir).glob("a_*"))

    # ensure only one of each file was found in the unzipped file
    assert 1 == len(i_file_glob) == len(s_file_glob), f"Error: there should be one and only 1 i_*, and s_* file"

    # pull the single file from the list of one
    i_file = i_file_glob[0]
    s_file = s_file_glob[0]

    # parse into tables (study and assay files) or dictionary of tables (investigation file)
    i_df_dict = isatab.load_investigation(fp=i_file.open())
    s_df = isatab.load_table(fp=s_file.open())
    a_dfs = [isatab.load_table(fp=a_file) for a_file in a_file_glob]


#def get_assay_type_from_i_df_assays()

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

def _unzip_ISA(isa_zip_path: str) -> str:
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

def parse_isa_dir_from_zip(isa_zip_path: str, pretty_print: bool = False) -> ISATabRecord:
    """ Unzips ISA zip files as found in GLDS metadata folders.

    ISA record metadata signature should match the specs here:
    https://isa-specs.readthedocs.io/en/latest/isamodel.html

    :param isa_zip_path: path to isa zip file
    :param pretty_print: print contents of parsed file, useful for debugging
    """
    INVESTIGATION_INDENT = 1
    STUDY_INDENT = 2
    ASSAY_INDENT = 3

    isa_temp_dir = _unzip_ISA(isa_zip_path)

    try:
        investigation = isatab_parser.parse(isatab_ref=isa_temp_dir)
    except UnicodeDecodeError as e:
        traceback.print_exc()
        raise Exception(f"ERROR: ISA.ZIP contents likely includes non-utf-8 decodable characters.  In the past, these includes things like degree signs in temperatures.")

    # only print if requested, useful for debugging parsing and data extraction
    if pretty_print:
        print(f"INVESTIGATION: ISA ZIP: {os.path.basename(isa_zip_path)}")
        print("="*95)
        _print_metadata(investigation, level=INVESTIGATION_INDENT)
        for i, study in enumerate(investigation.studies):
            print("\n")
            print(f"STUDY {i+1} of {len(investigation.studies)}")
            print("="*95)
            _print_metadata(study, level=STUDY_INDENT)
            for design_descriptor in study.design_descriptors:
                [print("\t"*STUDY_INDENT + f"{k}:  {v}") for k,v in design_descriptor.items()]

            # iterature thorugh assays
            for j, assay in enumerate(study.assays):
                print("\n")
                print(f"ASSAY {j+1} of {len(study.assays)} from STUDY {i+1}")
                print("="*95)
                _print_metadata(assay, level=ASSAY_INDENT)

    return investigation

class AssayNotFoundException(Exception):
    pass

def get_assay(study,
              assay_measurement_type: str,
              assay_technology_type: str):
    """ Returns an assay that matches the supplied measurement and technology type


    Matching based on regex
        - is case insensitive
        - treats the following between words as equivalent (" ","-","_")
        - No seperator between words is also accepted

    :param assay_measurement_type: Measurement type.  E.G.transcription profiling
    :param assay_technology_type: Technology type.  E.G. RNA Sequencing (RNA-Seq)
    """
    measurement_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_measurement_type).split("\ "))
    technology_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_technology_type).split("\ "))

    found_assays = list()

    for i, assay in enumerate(study.assays):
        found_assays = (i, assay.metadata['Study Assay Measurement Type'], assay.metadata['Study Assay Technology Type'])
        if re.match(measurement_regex, assay.metadata['Study Assay Measurement Type']) and\
           re.match(technology_regex, assay.metadata['Study Assay Technology Type']):
            return assay
    else:
        print(f"Found Assays: {found_assays}")
        raise AssayNotFoundException(f"Did not find compatible assay after scanning metadata for measurement type: {assay_measurement_type} and technology type: {assay_technology_type}")

def query_assay_from_a_dfs(a_dfs,
                        assay_measurement_type: str,
                        assay_technology_type: str):
    measurement_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_measurement_type).split("\ "))
    technology_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_technology_type).split("\ "))

    for i, assay_df in enumerate(a_dfs):
        print(assay_df.columns)
        found_measurement = assay_df['Study Assay Measurement Type'].unique()[0]
        found_tech_type = assay_df['Study Assay Technology Type'].unique()[0]
        found_assays = (i, found_measurement, found_tech_type)
        if re.match(measurement_regex, found_measurement) and\
           re.match(technology_regex, found_tech_type):
            return assay_df
    else:
        print(f"Found Assays: {found_assays}")
        raise AssayNotFoundException(f"Did not find compatible assay after scanning metadata for measurement type: {assay_measurement_type} and technology type: {assay_technology_type}")

def query_assay_from_i_df_assays(i_df_assays,
                                assay_measurement_type: str,
                                assay_technology_type: str):
    measurement_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_measurement_type).split("\ "))
    technology_regex = "(?i)\s*" + "[ ,\-,_]*".join(re.escape(assay_technology_type).split("\ "))

    mask_matches_measurement = i_df_assays['Study Assay Measurement Type'].str.match(measurement_regex)
    mask_matches_technology = i_df_assays['Study Assay Technology Type'].str.match(technology_regex)

    x_df = i_df_assays.loc[mask_matches_measurement & mask_matches_technology]
    if len(x_df) != 1:
        raise AssayNotFoundException(f"Did not find compatible assay after scanning metadata for measurement type: {assay_measurement_type} and technology type: {assay_technology_type}")
    else:
        # take this single row, transpose and return as dictionary
        return x_df.T.to_dict()[1]

def get_factor_names(study):
    return  [f"{study_factor['Study Factor Name']}"
            for study_factor in study.factors]

def extract_has_ercc(study):
    protocols = study.protocols
    spike_in_protocols = [protocol for protocol in protocols if protocol["Study Protocol Type"] == "spike-in quality control role"]
    if len(spike_in_protocols) == 1:
        return True
    elif len(spike_in_protocols) == 0:
        return False
    else:
        raise ValueError(f"Unexpectedly found more than 1 spike-in protocol. {spike_in_protocols}")

def extract_organism(study):
    """ returns organism studied.
    Assumption: each GLDS study is for one organism
    """
    try:
        for node, node_data in study.nodes.items():
            if "source" in node:
                if  organism_meta_attr := node_data.metadata.get("Characteristics[Organism]"):
                    if organism := getattr(organism_meta_attr[0], "Organism", None):
                        return organism
                # catch lowercase cases
                elif  organism_meta_attr := node_data.metadata.get("Characteristics[organism]"):
                    if organism := getattr(organism_meta_attr[0], "organism", None):
                        return organism
                # catch entries with mesh
                elif  organism_meta_attr := node_data.metadata.get("Characteristics[Organism,Http://Purl.Bioontology.Org/Ontology/Sty/T001,Mesh]]"):
                        if organism := getattr(organism_meta_attr[0], "organism", None):
                            return organism
    except Exception as e:
        print(e)
    raise ValueError(f"Could not find organism data. Last node metadata: {node_data.metadata}")

def extract_read_length(node_data: NodeRecord) -> tuple(int, int):
    """ Extract read length

    :param node_data: A node record to parse for read length metadata
    """
    read_lengths = None
    try:
        if  read_length_meta_attr := node_data.metadata.get("Parameter Value[Read Length,http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C153362,NCIT]"):
            if read_length := getattr(read_length_meta_attr[0], "Read_Length_http_ncicb_nci_nih_gov_xml_owl_EVS_Thesaurus_owl_C153362_NCIT", None):
                read_lengths =  (read_length, read_length)
        elif  read_length_meta_attr := node_data.metadata.get("Parameter Value[Read Length]"):
            if read_length := getattr(read_length_meta_attr[0], "Read_Length", None):
                read_lengths = (read_length, read_length)
        # catch studies with different read lengths
        elif  read_length_R1_meta_attr := node_data.metadata.get("Parameter Value[R1 Read Length]"):
            read_length_R2_meta_attr = node_data.metadata.get("Parameter Value[R2 Read Length]")
            if read_length_R1 := getattr(read_length_R1_meta_attr[0], "R1_Read_Length", None):
                read_length_R2 = getattr(read_length_R2_meta_attr[0], "R2_Read_Length", None)
                read_lengths = (read_length_R1, read_length_R2)

    except Exception as e:
        print(e)
    if read_lengths:
        return read_lengths
    else:
        raise ValueError(f"Could not find read length data. Last node metadata: {node_data.metadata}")

def get_glds_filelisting_json(accession: str) -> (dict, str):
    """ Return filelisting json accession number """
    # extract file urls
    glds_json = read_json(GLDS_URL_PREFIX + accession)
    try:
        _id = glds_json[0]["_id"]
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    file_list_url = FILELISTINGS_URL_PREFIX + _id
    file_listing_json = read_json(file_list_url)
    return file_listing_json, glds_json[0]["version"]

def get_organism_col_name(merged_df: DataFrame) -> str:
    def valid_organism_col(col_name: str):
        if any((
            col_name.startswith("Characteristics[Organism"),
            col_name.startswith("Characteristics[organism")
            )):
            return True
        return False

    organism_col = [col for col in merged_df.columns if valid_organism_col(col)]
    assert len(organism_col) == 1, f"Should be one and only one organism column, but found {len(organism_col)}: {organism_col}"
    return organism_col[0]

def clean_up_quotes(i_file: Path) -> Path:
    """ Converts double quotes within double quoted fields to single quotes.
    Also removes trailing tabs before newlines.
    """
    fixed_lines = list()
    cleaned_file_name = Path("cleaned_" + str(i_file.name))
    with open(cleaned_file_name, "w") as corrected_i_file:
        with open(i_file, "r") as f:
            for i, line in enumerate(f.readlines()):
                fixed_tokens = list()
                tokens = line.split('\t')
                for token in tokens:
                    if token == '\n':
                        # occurs when an additional trailing tab is included
                        continue
                    # assume double quoted field
                    if token[0] == '"':
                        token = '"' + token[1:-1].replace('"',"'") + '"'
                        fixed_lines.append(i+1)
                    fixed_tokens.append(token)
                corrected_line = '\t'.join(fixed_tokens)
                if corrected_line[-1] != '\n':
                    corrected_line = corrected_line + '\n'
                corrected_i_file.write(corrected_line)
    if fixed_lines:
        print("Fixed double quotes for lines: {}")
    return Path(cleaned_file_name)

def isa_to_RNASeq_runsheet(isazip, accession):
    isa = parse_isa_dir_from_zip(isazip)

    # extract study
    # assert only one study
    # ASSSUMPTION: GeneLab ISA files only contain one study
    assert len(isa.studies) == 1
    study = isa.studies[0]



    assay = get_assay(study,
                      assay_measurement_type = "transcription profiling",
                      assay_technology_type = "RNA Sequencing (RNA-Seq)")
    # start tracking project-wide data
    project = dict()
    project["GLDS"] = accession

    file_listing_json, project["version"] = get_glds_filelisting_json(accession)

    # extract samples from assay
    samples = dict()
    SAMPLE_PREFIX = "sample-"

    FILE_EXTRACTION_FUNCTIONS = {
        "Parameter Value[Merged Sequence Data File]":_extract_files_merged,
        "Raw Data File":_extract_file_raw,
                        }
    ALLOWED_RAW_READS_KEY = set(FILE_EXTRACTION_FUNCTIONS.keys())

    # extract factor value names from study
    factor_names = get_factor_names(study)

    # extract factor value values to samples
    for node_key, node_data in assay.nodes.items():
        if node_key.startswith(SAMPLE_PREFIX):
            sample_name = node_data.name.strip()
            print(f"Extracting data for sample: {sample_name}")
            samples[sample_name] = dict()
            samples[sample_name]["node"] = node_data
            study_node = study.nodes[node_key]

            samples[sample_name]["factors"] = dict()
            for factor_name in factor_names:
                #print(f"Factor: {factor_name}")
                factor_key = f"Factor Value[{factor_name}]"
                if node_data.metadata.get(factor_key):
                    # factor names as attributes replace spaces with '_'
                    parse_1 = getattr(node_data.metadata.get(factor_key)[0], factor_name.replace(" ","_").replace("-","_"), None)
                    # parse 2 example: GLDS-235 specifically key : Tissue Homogenate Preservation Time at -80C in RLT Buffer that maps to attribure Tissue_Homogenate_Preservation_Time_at_80C_in_RLT_Buffer
                    parse_2 = getattr(node_data.metadata.get(factor_key)[0], factor_name.replace(" ","_").replace("-","_").replace("__","_"), None)
                    factor_value = parse_1 if parse_1 else parse_2
                    #print(f"ASSAY: {factor_value}")
                    samples[sample_name]["factors"][factor_key] = factor_value
                if study_node.metadata.get(factor_key):
                    #print(study_node.metadata.get(factor_key)[0])
                    # factor names as attributes replace spaces with '_'
                    parse_1 = getattr(study_node.metadata.get(factor_key)[0], factor_name.replace(" ","_").replace("-","_"), None)
                    # parse 2 example: GLDS-235 specifically key : Tissue Homogenate Preservation Time at -80C in RLT Buffer that maps to attribure Tissue_Homogenate_Preservation_Time_at_80C_in_RLT_Buffer
                    parse_2 = getattr(study_node.metadata.get(factor_key)[0], factor_name.replace(" ","_").replace("-","_").replace("__","_"), None)
                    factor_value = parse_1 if parse_1 else parse_2
                    samples[sample_name]["factors"][factor_key] = factor_value
                #
                if factor_value == None:
                    print(node_key, node_data)
                    print("***")
                    print(study_node.metadata)
                    raise ValueError(f"A value MUST exist for each sample and each factor. {(sample_name, factor_key, factor_value)}")


            # extract filenames
            # find valid keys
            valid_keys = [key for key in node_data.metadata.keys() if key in ALLOWED_RAW_READS_KEY]
            if set(("Raw Data File", "Parameter Value[Merged Sequence Data File]")).issubset(set(valid_keys)):
                valid_keys = ["Parameter Value[Merged Sequence Data File]"]

            # only one valid key should be associated
            assert len(valid_keys) == 1, f"Found keys: {node_data.metadata.keys()}"
            isa_key = valid_keys[0]
            project['isa_key'] = isa_key
            file_names = FILE_EXTRACTION_FUNCTIONS[isa_key](node_data)

            # There should be 1 or 2 file names per sample only
            # check and set for each sample
            assert len(file_names) in (1,2), f"Unexpected number of file_names ({len(file_names)}) for {sample_name}. File names {file_names}"
            if len(file_names) == 1:
                project['paired_end'] = False
            elif len(file_names) == 2:
                project['paired_end'] = True
            samples[sample_name]["file_names"] = file_names

            # extract read length: tuple(R1, R2) or if single tuple(R1, R1)
            samples[sample_name]["read_length"] = extract_read_length(node_data)

            file_urls = list()
            for file_name in file_names:
                #print(file_name)
                # uses alternative url
                valid_urls = [f"{GENELAB_ROOT}/genelab/static/media/dataset/{quote(entry['file_name'])}?version={entry['version']}" for entry in file_listing_json if entry["file_name"] == file_name]

                # only one valid url should be associated
                assert len(valid_urls) == 1, f"Bad Number Of Listings For File: {file_name}. There are {len(valid_urls)} listings for this file. There should be one and only one entry. File listing URL: {file_list_url}"
                file_urls.append(valid_urls[0])
            samples[sample_name]["file_urls"] = file_urls

    # extract additional project wise metadata
    project["has_ercc"] = extract_has_ercc(study)
    project["organism"] = extract_organism(study)

    return write_proto_runsheet(accession, samples, project)

def isa_to_Microarray_runsheet(isazip, accession, missing_col_allowed=False):
    isa_temp_dir = _unzip_ISA(isazip)
    # find files using glob
    i_file_glob = list(Path(isa_temp_dir).glob("i_*"))
    s_file_glob = list(Path(isa_temp_dir).glob("s_*"))
    a_file_glob = list(Path(isa_temp_dir).glob("a_*"))

    # ensure only one of each file was found in the unzipped file
    assert 1 == len(i_file_glob) == len(s_file_glob), f"Error: there should be one and only 1 i_*, and s_* file"

    # pull the single file from the list of one
    i_file = i_file_glob[0]
    s_file = s_file_glob[0]

    # correct i_file
    i_file = clean_up_quotes(i_file)
    print(i_file)

    # parse into tables (study and assay files) or dictionary of tables (investigation file)
    i_df_dict = isatab.load_investigation(fp=i_file.open())
    s_df = isatab.load_table(fp=s_file.open())

    # isolate correct assay
    assert len(i_df_dict["s_assays"]) == 1, "Should be length of one (dataframe)"
    i_assay_dict = query_assay_from_i_df_assays(
                                    i_df_dict["s_assays"][0],
                                    assay_measurement_type = "transcription profiling",
                                    assay_technology_type = "DNA Microarray"
                                    )
    print(i_assay_dict)
    # Load assay table
    assay_file = isa_temp_dir / Path(i_assay_dict["Study Assay File Name"])
    print(f"Loading Assay File: {assay_file}")
    a_df = isatab.load_table(fp=assay_file.open())

    # merge a_df and s_df
    runsheet_df = a_df.merge(s_df, on="Sample Name", suffixes=("_assay","_sample"))
    assert len(runsheet_df) != 0, f"Empty runsheet after merge, check assay and study files for sample name parity"

    # titlecase all columns to make name matching easier
    runsheet_df.columns = map(str.title, runsheet_df.columns)

    # filter columns
    # i.e. drop columns that do not affect processing parameters
    KNOWN_DATA_FILE_COLUMN_VARIANTS = ["Array Data File","Parameter Value[array data file]"]
    array_data_column = [col for col in runsheet_df.columns if col in KNOWN_DATA_FILE_COLUMN_VARIANTS]
    assert len(array_data_column) == 1, f"One and only one column that indicates the array data file should exist"

    factor_columns = [col for col in runsheet_df.columns if col.startswith("Factor Value[")]

    # organism column name
    organism_column = get_organism_col_name(runsheet_df)

    columns_to_keep = ["Sample Name",
                       "Source Name",
                       "Label",
                       "Hybridization Assay Name",
                       organism_column] + factor_columns + array_data_column

    try:
        missing_columns = set(columns_to_keep).difference(set(runsheet_df.columns))
        assert len(missing_columns) == 0, "Missing Required Metadata Columns Detected"
        runsheet_df = runsheet_df[columns_to_keep]
    except AssertionError: # this indicates a required column was not found
        if missing_col_allowed:
            print(f"MISSING COLUMNS ALLOWED: Missing columns: {missing_columns}.  Columns found: {list(runsheet_df.columns)}")
            # Create missing columns
            for missing_col in missing_columns:
                runsheet_df[missing_col] = "COULD NOT BE EXTRACTED FROM ISA"
            runsheet_df = runsheet_df[columns_to_keep]
        else:
            raise KeyError(f"Missing columns Error: {missing_columns}.  Columns found: {list(runsheet_df.columns)}")
    # rename array data column
    runsheet_df = runsheet_df.rename(mapper={variant:"array_data_file" for variant in KNOWN_DATA_FILE_COLUMN_VARIANTS}, axis='columns')

    # populate with additional dataset-wide
    runsheet_df["GLDS"] = accession
    runsheet_df["Study Assay Measurement Type"] = i_assay_dict["Study Assay Measurement Type"]
    runsheet_df["Study Assay Technolgy Type"] = i_assay_dict["Study Assay Technology Type"]
    runsheet_df["Study Assay Technology Platform"] = i_assay_dict["Study Assay Technology Platform"]

    # get file urls
    file_listing_json, runsheet_df["version"] = get_glds_filelisting_json(accession)


    def _get_file_url(file_name, file_listing_json):
        file_url = [f"{GENELAB_ROOT}/genelab/static/media/dataset/{quote(entry['file_name'])}?version={entry['version']}" for entry in file_listing_json if entry["file_name"] == file_name]
        assert len(file_url) == 1, f"One and only one file url should exist for each file. {file_name}"
        return file_url[0]
    # populate file urls
    runsheet_df["array_data_file_path"] = runsheet_df["array_data_file"].apply(lambda file_name: _get_file_url(file_name, file_listing_json))
    runsheet_df["is_array_data_file_compressed"] = runsheet_df["array_data_file"].str.endswith(".gz")

    tmp_proto = f"{accession}_proto_run_sheet.csv"
    runsheet_df.to_csv(tmp_proto, index=False)
    return tmp_proto

def write_proto_runsheet(accession: str, samples: dict, project: dict):
    output_file = Path(f"{accession}_RNASeq_runsheet.csv")
    with open(output_file, "w") as f:
        ###########################################################
        # WRITE HEADER ############################################
        ###########################################################
        factor_string = ','.join(samples[list(samples)[0]]['factors'].keys())
        f.write(f"sample_name,read1_url,"\
                f"paired_end,has_ERCC,version,organism,read_length_R1,isa_key,"\
                f"protocol,raw_read1,trimmed_read1,STAR_Alignment,RSEM_Counts,"\
                f"raw_read_fastQC,trimmed_read_fastQC,{factor_string}")
        if project["paired_end"]:
            f.write(",read2_url,raw_read2,trimmed_read2,read_length_R2")
        f.write("\n")

        ###########################################################
        # WRITE SAMPLE ROWS #######################################
        ###########################################################
        for sample_name, sample in samples.items():
            sample_factor_values_string = ','.join(sample['factors'].values())
            if project['paired_end']:
                read1 = sample["file_names"][0]
                read2 = sample["file_names"][1]
                read1_url = sample["file_urls"][0]
                read2_url = sample["file_urls"][1]
                read1_read_length = sample["read_length"][0]
                read2_read_length = sample["read_length"][1]
            else:
                read1 = sample["file_names"][0]
                read2 = None
                read1_url = sample["file_urls"][0]
                read2_url = None
                read1_read_length = sample["read_length"][0]
                read2_read_length = None

            f.write(f"{sample_name.replace(' ','_')},{read1_url},"\
                    f"{project['paired_end']},{project['has_ercc']},"\
                    f"{project['version']},{project['organism']},"\
                    f"{read1_read_length}"\
                    f",{project['isa_key']},anySampleType,raw_read1,"\
                    f"trimmed_read1,STAR_Alignment,RSEM_Counts,"\
                    f"raw_read_fastQC,trimmed_read_fastQC,"\
                    f"{sample_factor_values_string}")
            if project["paired_end"]:
                f.write(f",{read2_url},raw_read2,"\
                        f"trimmed_read2,{read2_read_length}")
            f.write("\n")

    return output_file

def _extract_files_merged(node):
    #print("Merged Extracted")
    files_string =  node.metadata["Parameter Value[Merged Sequence Data File]"][0].Merged_Sequence_Data_File
    return [file.strip() for file in files_string.split(",")]

def _extract_file_raw(node):
    #print("Raw Data File Extraction")
    files_string =  node.metadata["Raw Data File"][0]
    return [file.strip() for file in files_string.split(",")]

def main():
    args = _parse_args()
    isazip = download_isa(args.accession, args.alternate_url)

    #if args.to_runsheet:
    #    assay

    if args.to_RNASeq_runsheet:
        # generate proto run sheet from ISA
        proto_run_sheet = isa_to_RNASeq_runsheet(isazip, args.accession)
        shutil.copy(proto_run_sheet, "tmp_proto_run_sheet.csv")
        # load peppy project config
        if args.alt_peppy_template:
            with importlib.resources.path("AST", f"RNASeq_RCP_alt{args.alt_peppy_template}.yaml") as template:
                template_path = template
        else:
            with importlib.resources.path("AST", "RNASeq_RCP.yaml") as template:
                template_path = template
        shutil.copy(template_path, ".")
        p = peppy.Project(template_path.name)
        filled_run_sheet_name = f"AST_autogen_template_{template_path.with_suffix('').name}_{proto_run_sheet}"
        p.sample_table.to_csv(filled_run_sheet_name)
        print(f"Autogenerating Paths for RNASeq run sheet")
        print(f"Template (in AST package): {template_path.name}")
        print(f"Filled Run Sheet: {filled_run_sheet_name}")
        os.remove(proto_run_sheet)
        os.remove("tmp_proto_run_sheet.csv")

    elif args.to_Microarray_runsheet:
        # generate proto run sheet from ISA
        proto_run_sheet = isa_to_Microarray_runsheet(isazip, args.accession, args.allow_missing_columns)
        '''shutil.copy(proto_run_sheet, "tmp_proto_run_sheet.csv")
        # load peppy project config
        with importlib.resources.path("AST", "Microarray.yaml") as template:
            template_path = template
        shutil.copy(template_path, ".")
        p = peppy.Project(template_path.name)
        filled_run_sheet_name = f"AST_MICROARRAY_v1_{proto_run_sheet}"
        p.sample_table.to_csv(filled_run_sheet_name)
        print(f"Autogenerating Paths for Microarray run sheet")
        print(f"Template (in AST package): {template_path.name}")
        print(f"Filled Run Sheet: {filled_run_sheet_name}")
        os.remove(proto_run_sheet)
        os.remove("tmp_proto_run_sheet.csv")'''

    else:
        print(f"No Runsheet generation requested")

if __name__ == "__main__":
    main()
