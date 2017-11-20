#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.version import __version__

from setuptools import setup, find_packages


url = 'https://github.com/HadrienG/InSilicoSeq'

setup(
    name='InSilicoSeq',
    version=__version__,

    description='a sequencing simulator',

    url=url,
    download_url=url + '/tarball/' + __version__,
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
