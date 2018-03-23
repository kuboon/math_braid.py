from setuptools import setup, find_packages

setup(
	name = 'math_braid',
	packages = find_packages(exclude=['tests', 'tests.*']),
	version = '0.1',
	description = 'Pure python braid group implementation based on sympy.',
	author = 'kuboon',
	author_email = 'kuboon@trick-with.net',
	url = 'https://github.com/kuboon/braid_group.py',
	download_url = 'https://github.com/kuboon/braid_group.py',
	keywords = ['math', 'braid', 'group'],
	classifiers = [],
)
