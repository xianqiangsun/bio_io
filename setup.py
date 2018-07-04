#!/usr/bin/env python

from distutils.core import setup
from setuptools import setup, Extension, find_packages
import os

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(name='BioIO',
      version='1.0.1',
      description='Python Distribution Utilities',
      author='Xianqiang Sun',
      author_email='xianqiangleo@gmail.com',
      url='github.com/xianqiangsun/bio_io',
      packages=find_packages(),
      package_dir={'bio_io': './BioIO'},
      #package_data={
      #            'Bio': ['database/*']
      #            },
      #packages=['quantitative'],

      include_package_data=True,
      install_requires=['numpy',
                        'pandas',
                        'scipy',
                        'matplotlib',
                        'bioservices'
                        ])
