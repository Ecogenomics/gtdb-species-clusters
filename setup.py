from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'genome_tree_tk', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='GenomeTreeTK',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['genome_tree_tk'],
    scripts=['bin/genome_tree_tk'],
    package_data={'genome_tree_tk' : ['VERSION'], '': ['distributions/*.txt']},
    url='http://pypi.python.org/pypi/genome_tree_tk/',
    license='GPL3',
    description='A toolbox for working with genome trees.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.8.0",
        "biolib >= 0.0.1"],
)
