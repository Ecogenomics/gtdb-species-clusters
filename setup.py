from distutils.core import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'genometreetk', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='GenomeTreeTk',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['genometreetk', 'genometreetk.markers'],
    scripts=['bin/genometreetk'],
    package_data={'genometreetk' : ['VERSION'], '': ['distributions/*.txt']},
    url='http://pypi.python.org/pypi/genometreetk/',
    license='GPL3',
    description='A toolbox for working with genome trees.',
    install_requires=[
        "numpy >= 1.8.0",
        "biolib >= 0.0.45",
        "dendropy >= 4.0.0"],
)
