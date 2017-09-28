#!/usr/bin/python3.5
"""
Arthur Zwaenepoel

Version 1.1:
    - Rewrite of Ks distribution construction
    - Implementation of one-vs-one ortholog Ks distributions
"""

from setuptools import setup

setup(
    name='wgd',
    version='1.1',
    packages=['wgd'],
    url='http://github.com/arzwa/wgd',
    license='',
    author='Arthur Zwanepoel',
    author_email='arzwa@psb.vib-ugent.be',
    description='MORPH bulk CLI',
    py_modules=['wgd_cli'],
    include_package_data=True,
    install_requires=[
        'click',
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
    ],
    entry_points='''
        [console_scripts]
        wgd=wgd_cli:cli
    ''',
)
