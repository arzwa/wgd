#!/usr/bin/python3.5
"""
Arthur Zwaenepoel

Version 1.1:
    - Rewrite of Ks distribution construction
    - Implementation of one-vs-one ortholog Ks distributions
    - Phylogenetic trees for weighting (FastTree & PhyML)

Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be
"""

from setuptools import setup

setup(
    name='wgd',
    version='1.1',
    packages=['wgd'],
    url='http://github.com/arzwa/wgd',
    license='GPL',
    author='Arthur Zwanepoel',
    author_email='arzwa@psb.vib-ugent.be',
    description='MORPH bulk CLI',
    py_modules=['wgd_cli'],
    include_package_data=True,
    install_requires=[
        'click',
        'seaborn',
        'coloredlogs',
        'fastcluster',
        'peakutils',
        'numpy',
        'sklearn',
        'scipy',
        'matplotlib',
        'plumbum',
        'fastcluster',
        'pandas',
        'progressbar2',
        'joblib',
        'ete3',
    ],
    entry_points='''
        [console_scripts]
        wgd=wgd_cli:cli
    ''',
)
