from setuptools import setup, find_packages

setup(
    name='lisflood-preprocessing',
    version='0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'lfcoordinates=lisfloodpreprocessing.correct_coordinates:main',
        ],
    },
    install_requires=[
        'geopandas',
        'numpy',
        'pandas',
        'tqdm',
        'pyflwdir',
        'pyyaml',
        'rioxarray',
        'xarray',
    ],
    author='Peter Burek, Jesús Casado Rodríguez',
    author_email='burek@iiasa.ac.at, jesus.casado-rodriguez@ec.europa.eu',
    description='Package to preprocess inputs of the hydrological model LISFLOOD.',
    keywords='hydrology lisflood stations',
)