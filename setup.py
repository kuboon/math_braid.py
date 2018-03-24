from setuptools import setup, find_packages

setup(
	name='math_braid',
	packages=["math_braid", "math_braid.extras"],
	version='0.3',
	description='Pure python braid group implementation based on sympy.',
	long_description='''Pure python braid group implementation based on sympy.
	Original code is here: http://www-math.mit.edu/~seidel/geng/
	''',
	license='BSD',
	author='kuboon',
	author_email='kuboon@trick-with.net',
	url='https://github.com/kuboon/math_braid.py',
	keywords="Math Braid Permutation",
	python_requires='>=3',
	install_requires=["sympy>=1.0"],
	classifiers=[
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
)
