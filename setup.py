from setuptools import setup
#import warnings
#from pkg_resources import DeprecationWarning as pkg_warning
#warnings.simplefilter("ignore", pkg_warning)

setup(
    name='bio-allo',
    version='1.2.0',
    author='Alexis Morrissey',
    author_email='anm5579@psu.edu',
    packages=['Allo'],
    python_requires='>=3.10',
    scripts=['Allo/allo'],
    url='https://github.com/seqcode/allo/archive/refs/tags/1.2.0.tar.gz',
    license='LICENSE.txt',
    description='A multi-mapped read rescue strategy for gene regulatory data',
    include_package_data=True,
    install_requires=[
        'numpy',
        'joblib',
        'tensorflow>=2.17',
        'pysam'
    ]
)

