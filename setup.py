from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='AmorphSim',
      version='0.01',
      description='Tool for simulating amorphous Materials',
      long_description=readme(),
      keywords='Simulations STEM Electron Microscopy Glass',
      url='https://github.com/CSSFrancis/AmorphSim',
      author='CSSFrancis',
      author_email='csfrancis@wisc.edu',
      liscense='MIT',
      packages=['AmorphSim',
                'AmorphSim.utils',],
      install_requires=['hyperspy >=1.5',
                        'numpy>=1.10,!=1.70.0',
                        'matplotlib',
                        'scipy',
                        'skimage'],
      zip_safe=False)