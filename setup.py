#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:00:26 2024

@author: yavuzb2
"""

import setuptools

# Read the contents of the README file
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CombinationTherapyTargets",  
    version="0.1.0",  
    author="Bengi Ruken Yavuz",
    author_email="bengi.yavuz@nih.gov",
    description="Protein targets for Combination Therapies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bengiruken/CombinationTherapyTargets/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
        install_requires=[
                'pandas==2.0.3',
                'scipy==1.11.1',
                'numpy==1.24.3',
                'gseapy==1.1.2',
                'matplotlib-inline==0.1.6',
                'seaborn==0.12.2',
                'networkx==3.2.1',
                ]
    )
    
    
    
    
    
    
    
  