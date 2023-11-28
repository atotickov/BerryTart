__author__ = 'atotickov'

from pathlib import Path
from os.path import join, dirname
from setuptools import setup, find_packages

dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'biopython',
                'lxml', 'beautifulsoup4', 'requests', 'gzip']

setup(name='BerryTart',
      version='0.1',
      packages=find_packages(),
      author='Azamat Totikov',
      author_email='a.totickov1@gmail.com',
      install_requires=dependencies,
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=list(map(str, sorted(Path('python_scr/').rglob("*.py")))))
