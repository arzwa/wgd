#!/usr/bin/python3.8
# contact: arzwa@psb.vib-ugent.be
# https://stackoverflow.com/questions/43658870/requirements-txt-vs-setup-py

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='wgd',
    version='2.0.0',
    packages=['wgd'],
    url='http://github.com/arzwa/wgd',
    license='GPL',
    author='Arthur Zwanepoel',
    author_email='arzwa@psb.vib-ugent.be',
    description='wgd',
    py_modules=['cli'],
    include_package_data=True,
    install_requires=[
        # fill in
    ],
    entry_points='''
        [console_scripts]
        wgd=cli:cli
    ''',
)
