""" Functions for 'cleaning' isazip files to ensure parsing via isatools works
"""
from pathlib import Path
import re

def parse_i_file(i_file_path: Path) -> dict:
    """ A function to parse i_file sections for correction
    """
    with i_file_path.open() as f:
        for line in f.readlines():
            yield line

def i_file_clean_field_characters(line):
    """ A function that removes replaces double quoted terms with
    """
    ...

def cleanup_double_quotes(file_path: Path) -> Path:
    """ Replaces double quoted fields using U+0022 with analogous forward and reverse quotes
    “ (U+201C) and ”(U+201D).  Ignores any { "\t" } which are likely to fall between fields in
    when double quote protected fields are used in the file.

    """
    def replace_func(match):
        string = match.group(0)
        core = string[1:-1]
        replacement = '“' + core + '”'
        print(f"Replacing {string} with {replacement}")
        return replacement
    replace_this = re.compile('(?<!^)(?<!\t)"[^\t]*"(?!\t)(?!$)')
    #replace_this = re.compile("Investigation")
    #replace_this = re.compile(r'Ground')
    with file_path.open() as f:
        for line in f.readlines():
            #print(line)
            print(replace_this.sub(replace_func, line))
    return "NOMATCH"
