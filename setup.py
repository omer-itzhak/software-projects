from distutils.command import install
from setuptools import setup, find_packages, Extension

setup(
    name = "spkmeans.py",
    author = "Sharon Tirza and Omer Itzhak",
    packages = find_packages(where='.'),
    install_requires = ['invoke'],
    classifiers = ["Programming Language :: Python :: Implementation ::CPython"],
    ext_modules = [Extension("mykmeanssp", sources=['spkmeans.c', 'spkmeansmodule.c'])]
)