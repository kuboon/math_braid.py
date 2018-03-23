from setuptools import setup, find_packages

requires = ["sympy>=1.0"]

setup(
	name = 'math_braid',
	packages = ["math_braid", "math_braid.extras"],
	version = '0.1',
	description = 'Pure python braid group implementation based on sympy.',
	author = 'kuboon',
	author_email = 'kuboon@trick-with.net',
	url = 'https://github.com/kuboon/braid_group.py',
	keywords = ['math', 'braid', 'group'],
	classifiers = [],
	python_requires='>=3',
)
