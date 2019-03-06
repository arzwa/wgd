#!/usr/bin/python3.5
"""
Arthur Zwaenepoel

Version 1.0

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

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='wgd',
    version='1.0',
    packages=['wgd'],
    url='http://github.com/arzwa/wgd',
    license='GPL',
    author='Arthur Zwanepoel',
    author_email='arzwa@psb.vib-ugent.be',
    description='wgd',
    py_modules=['wgd_cli'],
    include_package_data=True,
    install_requires=[
        'click>=7.0',
        'seaborn>=0.9.0',
        'coloredlogs>=10.0',
        'fastcluster==1.1.25',
        'numpy>=1.16',
        'sklearn',
        'scipy>=1.2',
        'matplotlib>=3.0.2',
        'plumbum>=1.6.7',
        'pandas==0.24.1',
        'progressbar2>=3.39',
        'joblib==0.11', # 0.12 seems to break the logging in parallel for loops
        'ete3>=3.1',
        'bokeh>=1.0.4'
    ],
    entry_points='''
        [console_scripts]
        wgd=wgd_cli:cli
    ''',
)
