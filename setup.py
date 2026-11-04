"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import os, re

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

VERSION_FILE = os.path.join(here, "rosetta", "version.py")

def get_version():
    with open(VERSION_FILE, "r") as f:
        content = f.read()
    match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', content, re.M)
    if not match:
        raise RuntimeError(f"Unable to find version string in {VERSION_FILE}")
    return match.group(1)

setup(
    name='scipion-chem-rosetta',  # Required
    version=get_version(),  # Required
    description='Scipion plugin in order to use tools provided by Rosetta software suite.',  # Required
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
