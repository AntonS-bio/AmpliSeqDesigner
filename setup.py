#!/usr/bin/env python3
from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


__version__ = 'v0.1_beta.10'

setup(name='HandyAmpliconTool',
      version=__version__,
      description='HandyAmpliconTool',
      long_description=readme(),
      python_requires='>=3.10.0',
      classifiers=['Development Status :: Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research'],
      keywords='PCR amplicon primers',
      url='https://github.com/AntonS-bio/HandyAmpliconTool.git',
      author='Anton Spadar',
      author_email='',
      packages=['scripts'],
      scripts=[
          'scripts/run.py'
      ]
)