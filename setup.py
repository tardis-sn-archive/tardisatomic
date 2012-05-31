#setup script for tardissatomic

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import glob

#version number
version = '0.01dev'

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
try:
	scripts.remove('scripts/README.rst')
except ValueError:
	pass


setup(name='tardisatomic',
    description='Atomic data for the TARDIS SN synthesis program',
    author='Wolfgang Kerzendorf (wkerzendorf@gmail.com)',
    version=version,
    packages=['tardisatomic'],
    package_data={'tardisatomic': ['data/*']},
    scripts=scripts
      )
      
