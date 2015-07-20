#!/usr/bin/env python

from setuptools import setup

setup(
    name='svfilter',
    version='0.0.1',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['svfilter = src.main:main']
    },
    url='https://github.com/bjpop/svfilter',
    license='LICENSE.txt',
    description='sv_filter is a program for filtering structural variants\
     overlapping with a set of genomic coordindates.',
    long_description=open('README.md').read(),
    install_requires=[
        "PyVCF == 0.6.7",
        "bx-python == 0.7.3",
    ],
)
