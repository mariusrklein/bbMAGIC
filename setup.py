#!/usr/bin/env python3

import setuptools

setuptools.setup(
    name="bbMAGIC",
    version="0.0.1",
    author="Marius Klein",
    author_email="marius.klein@embl.de",
    packages=setuptools.find_packages(),
    description="This package provides functions for the batch integration of single-cell omics data.",
    url="https://git.embl.de/mklein/bbmagic",
    license='MIT',
    python_requires='>=3.0',
    install_requires=[
        "scanpy>=1.9.1",
        "bbknn",
        "magic"
    ]
)