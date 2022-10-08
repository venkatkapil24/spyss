from setuptools import setup, find_packages

DESCRIPTION =  "Stochastic Pythonic Structure Search (SPySS)"
LONG_DESCRIPTION = 'README'

setup(
        name="spyss", 
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["numpy", "ase"], 
        ],
)
