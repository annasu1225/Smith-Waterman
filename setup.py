from setuptools import setup, find_packages

setup(
    name='SW',
    version='1.0.0',
    url='https://github.com/annasu1225/Smith-Waterman',
    author='Anna Su',
    author_email='anna.su@yale.edu',
    description='Smith-Waterman local alignment using blosum62 with open and extend gap penalties',
    packages=find_packages(),    
)
