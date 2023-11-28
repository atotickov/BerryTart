#!/usr/bin/env python3

from pathlib import Path
from os.path import join, dirname
from setuptools import setup, find_packages

dependencies = ['argparse', 'collections', 'gzip']

setup(name='BerryTart',
      packages=find_packages(),
      install_requires=dependencies,
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=list(map(str, sorted(Path('python_scr/').rglob("*.py")))))
