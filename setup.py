from setuptools import setup

setup(
	name='Allo',
	version='1.0',
	author='Alexis Morrissey',
	author_email='anm5579@psu.edu',
	packages=['Allo'],
	python_requires='>=3.10, <3.11',
	scripts=['Allo/allo'],
	url='https://github.com/seqcode/allo/',
	license='LICENSE.txt',
	description='A multi-mapped read rescue strategy for ChIP-seq data',
	install_requires=[
		'numpy >= 1.23.1',
    'joblib >= 1.1.0',
    'tensorflow == 2.11',
    'pysam == 0.20.0'
	]
)