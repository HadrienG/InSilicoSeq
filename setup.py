#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup

from iss.version import __version__

url="https://github.com/stavperez/InSilicoSeq"

with open("README.md") as f:
    long_description = f.read()

setup(
    name="InSilicoSeq",
    version=__version__,
    description="an annotation-aware fork of the InSilicoSeq sequencing simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=url,
    download_url=url + "/tarball/" + __version__,
    author="Hadrien Gourlé",
    author_email="hadrien.gourle@slu.se",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.10",
    tests_require=["pytest"],
    install_requires=["numpy", "scipy", "biopython>=1.80", "pysam>=0.15.1", "requests", "pandas", "gffutils"],
    include_package_data=True,
    entry_points={
        "console_scripts": ["iss = iss.app:main"],
    },
)
