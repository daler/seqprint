import ez_setup
ez_setup.use_setuptools()

import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'seqprint', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')

long_description = """
Extendible pretty-printing of genomic sequences, with support for various ASCII
annotations
"""

setup(
        name="seqprint",
        version=version,
        install_requires=['pybedtools', 'biopython', 'numpy'],
        packages=['seqprint',
                  'seqprint.test',
                  'seqprint.test.data'],
        author="Ryan Dale",
        description='pretty-printing of genomic sequences',
        long_description=long_description,
        url="none",
        package_data = {'seqprint':["test/data/*"]},
        package_dir = {"seqprint": "seqprint"},
        scripts = [],
        author_email="dalerr@niddk.nih.gov",
        classifiers=['Development Status :: 4 - Beta'],
    )
