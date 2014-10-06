#!/usr/bin/env python

from distutils.core import setup

setup(
    name='hitime',
    version='1.0.0',
    author='Andrew Isaac',
    author_email='aisaac@unimelb.edu.au',
    packages=['hitime'],
    entry_points={
        'console_scripts': ['hitime = hitime.hitime:main']
    },
    url='https://github.com/bjpop/HiTIME',
    license='LICENSE.txt',
    description=(
        'HiTIME: High-resolution Twin-Ion Metabolite Extraction.'),
    long_description=(
        'HiTIME searches for twin-ion pairs in high resolution Liquid Chromatography '
        'Mass Spectrometry (LCMS) data.'),
    install_requires=[
        "numpy == 1.7.1",
        "lxml == 3.2.3",
        "pymzml == 0.7.4",
    ],
)
