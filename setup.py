from setuptools import setup, find_packages

def read_version():
    with open("VERSION") as version_file:
        return version_file.read().strip()
    
setup(
    name='lisflood-preprocessing',
    version=read_version(),
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'lfcoords=lisfloodpreprocessing.main:main',
        ],
    },
    install_requires=[
        'geopandas',
        'numpy',
        'pandas',
        'tqdm',
        'pyflwdir',
        'shapely>=2.0',
        'pyyaml',
        'rioxarray',
        'xarray',
    ],
    author='Peter Burek, Jesús Casado Rodríguez',
    author_email='burek@iiasa.ac.at, chus.casado.88@gmail.com',
    description='Package to preprocess inputs of the hydrological model LISFLOOD.',
    keywords='hydrology lisflood stations',
)