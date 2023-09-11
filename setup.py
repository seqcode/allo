from setuptools import setup

setup(
	name='bio-allo',
	version='1.0.5',
	author='Alexis Morrissey',
	author_email='anm5579@psu.edu',
	packages=['Allo'],
	python_requires='>=3.0',
	scripts=['Allo/allo'],
	url='https://github.com/seqcode/allo/archive/refs/tags/v1.0.5.tar.gz',
	license='LICENSE.txt',
	description='A multi-mapped read rescue strategy for ChIP-seq data',
	include_package_data = True,
        install_requires=[
	'numpy',
        'joblib',
        'tensorflow >= 2.11',
        'pysam'
	]
)
