#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="selfies",
    version="1.0.2",
    author="Mario Krenn",
    author_email="mario.krenn@utoronto.ca, alan@aspuru.com",
    description="SELFIES (SELF-referencIng Embedded Strings) is a "
                "general-purpose, sequence-based, robust representation of "
                "semantically constrained graphs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aspuru-guzik-group/selfies",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5'
)
