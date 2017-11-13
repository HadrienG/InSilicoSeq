#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='InSilicoSeq',
    version='0.8.1',

    description='a sequencing simulator',

    url='https://github.com/HadrienG/InSilicoSeq',
    download_url='https://github.com/HadrienG/InSilicoSeq/tarball/0.8.1',
    author='Hadrien Gourlé',
    author_email='hadrien.gourle@slu.se',

    license='MIT',
    packages=find_packages(),

    tests_require=['nose'],
    install_requires=['numpy', 'scipy', 'biopython', 'pysam', 'future',
                      'joblib'],
    include_package_data=True,

    entry_points={
        'console_scripts': ['iss = iss.app:main'],
    }
)
