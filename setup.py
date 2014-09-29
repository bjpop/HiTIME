#!/usr/bin/env python

from distutils.core import setup

setup(
    name='HiTIME',
    version='1.0.0',
    author='Andrew Isaac',
    author_email='aisaac@unimelb.edu.au',
    packages=['HiTIME'],
    scripts=[''],
    entry_points={
        'console_scripts': ['hitme = hitime.hitime:main']
    },
    url='https://github.com/bjpop/HiTIME',
    license='LICENSE.txt',
    description=(
        'HiTIME: High-resolution Twin-Ion Metabolite Extraction.'),
    long_description=(
        'HiTIME searches for twin-ion pairs in high resolution Liquid Chromatography '
        'Mass Spectrometry (LCMS) data.'),
    install_requires=[
        "numpy >= 1.7.1",
        "scipy >= 0.12.0",
        "lxml >= 3.2.3",
        "pymzml >= 0.7.4",
        "Rtree >= 0.7.0",
        "matplotlib >= 1.3.1"
    ],
)
