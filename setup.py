#!/usr/bin/env python 

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name    = "selfies",
    version = "0.1.2",
    author  = "Mario Krenn",
    author_email = "mario.krenn@utoronto.ca, alan@aspuru.com",
    description  = "Self-referencing embedded strings",
    long_description = long_description,
    long_description_content_type = "text/markdown",
	url = "https://github.com/aspuru-guzik-group/selfies",
    packages = setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)

