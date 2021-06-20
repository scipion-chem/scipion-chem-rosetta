"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-rosetta',  # Required
    version='0.1',  # Required
    description='Scipion plugin in order to use tools provided by Rosetta software suite. Currently it has protocols for the use of the Rosetta DARC docking software',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-chem/scipion-chem-rosetta',  # Optional
    author='Alberto M. Parra-Perez',  # Optional
    author_email='amparraperez@gmail.com',  # Optional
    keywords='scipion docking ppi scipion-3.0 cheminformatics',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'rosetta = rosetta'},
    package_data={  # Optional
       'rosetta': ['rosetta_icon.png', 'protocols.conf'],
    }
)
