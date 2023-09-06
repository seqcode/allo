from setuptools import setup

setup(
	name='bio-allo',
	version='1.0.2',
	author='Alexis Morrissey',
	author_email='anm5579@psu.edu',
	packages=['Allo'],
	python_requires='>=3.10, <3.11',
	scripts=['Allo/allo'],
	url='https://github.com/seqcode/allo/archive/refs/tags/v1.0.2.tar.gz',
	license='LICENSE.txt',
	description='A multi-mapped read rescue strategy for ChIP-seq data',
	include_package_data = True,
        install_requires=[
	'numpy == 1.23',
        'joblib >= 1.1.0',
        'tensorflow == 2.11',
        'pysam == 0.20.0'
	]
)
