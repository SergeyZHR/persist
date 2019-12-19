#!/usr/bin/python3
# -*- coding: UTF-8 -*-

from setuptools import setup, find_packages
from os.path import join, dirname
import persist
setup(
	name='persist',
	version=persist.__version__,
	packages=find_packages(),
	entry_points={
		'console_scripts':['persist = persist.persist:main']
	},
	include_package_data=True,
	test_suite='tests',
	install_requires=['numpy', 'astropy', 'pandas', 'argparse'],
	long_description = open(join(dirname(__file__), 'README.md')).read(),
)