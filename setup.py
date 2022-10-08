from setuptools import setup, find_packages

DESCRIPTION =  "Stochastic Pythonic Structure Search (SPySS)"
LONG_DESCRIPTION = 'README'

setup(
        name="spyss", 
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
install_requires = [
    'ase',
    'numpy>=1.17.0',  # July 2019
    'scipy>=1.3.1',  # August 2019
    'matplotlib>=3.1.0',  # May 2019
    'importlib-metadata>=0.12;python_version<"3.8"'
]
)
