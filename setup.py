import os
from setuptools import setup


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setupDir, 'gtdb_species_clusters', 'VERSION')) as versionFile:
        return versionFile.readline().strip()


setup(
    name='gtdb_species_clusters',
    version=version(),
    python_requires='>=3.7',
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['gtdb_species_clusters'],
    scripts=['bin/gtdb_species_clusters'],
    package_data={'gtdb_species_clusters': ['VERSION']},
    url='https://github.com/Ecogenomics/gtdb-species-clusters',
    license='GPL3',
    description='This toolkit provides functionality for establishing, updating, and validating the species clusters used in the Genome Taxonomy Database.',
    install_requires=[
        "numpy >= 1.8.0",
        "dendropy >= 4.0.0",
        "gtdblib >= 1.10.1",
        "Levenshtein",
        "tqdm",
        "toml; python_version < '3.11'"],
)
