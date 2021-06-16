import codecs
import os
from setuptools import setup

dirname = os.path.dirname(__file__)

setup(
   name='AST',
   version='0.3.0',
   description='Tools for staging and setting up samplesheets.',
   author='Jonathan Oribello',
   author_email='jonathan.d.oribello@gmail.com',
   packages=['AST'],  #same as name
   scripts=[
            'AST/retrieve_isa_from_genelab.py',
           ],
   include_package_data=True,
   package_data={
        "": ['RNASeq_RCP.yaml', 'RNASeq_RCP_alt1.yaml', 'Microarray.yaml'],
                },
   setup_requires=['pytest-runner'],
   tests_require=['pytest']
)
