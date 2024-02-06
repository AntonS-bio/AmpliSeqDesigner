#!/usr/bin/env python3
from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
#exec(open('kleborate/version.py').read())

setup(name='HandyAmpliconTool',
      version=__version__,
      description='HandyAmpliconTool',
      long_description=readme(),
      classifiers=['Development Status :: Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research'],
      keywords='PCR amplicon primers',
      url='https://github.com/AntonS-bio/HandyAmpliconTool.git',
      author='Anton Spadar',
      author_email='',
      packages=['HandyAmpliconTool'],
      install_requires=['biopython>=1.81', 'ncbi-blast>=2.14.1'],
      scripts=[
          'scripts/run.py'
      ]
)