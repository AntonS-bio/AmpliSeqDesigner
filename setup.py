#!/usr/bin/env python3
import os
import sys
import subprocess
from setuptools import setup
from setuptools.command.install import install


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
#exec(open('kleborate/version.py').read())


def check_dir_write_permission(directory):
    if os.path.isdir(directory) and not os.access(directory, os.W_OK):
        sys.exit('Error: no write permission for ' + directory + '  ' +
                 'Perhaps you need to use sudo?')


def build_blast_db(data_dir, fasta_filename, seq_type):
    fasta_path = os.path.join(data_dir, fasta_filename)
    makeblastdb_cmd = ['makeblastdb', '-dbtype', seq_type, '-in', fasta_path]
    print('  ' + ' '.join(makeblastdb_cmd))
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(makeblastdb_cmd, stdout=devnull)


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
      install_requires=['biopython', 'ncbi-blast'],
      entry_points={'console_scripts': ['handyamplicontool = run:main']},
      include_package_data=True,
      zip_safe=False,
      cmdclass={'install': none})