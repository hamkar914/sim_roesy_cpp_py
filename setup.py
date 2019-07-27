from setuptools import setup, Extension
import numpy as np

module = Extension("NMR_utils",sources=["NmrModule.cpp"], language='c++',extra_compile_args=["-std=c++11", "-Wall", "-O3"], include_dirs=["inc/"])

setup(
      name="NMR_utils",
      ext_modules = [module],
      include_dirs = [np.get_include()],
      version="1.0",
      author="hamkar",
      author_email="hkarlsson914@gmail.com",
      description="Simulate cross relaxation during ROESY/NOESY mixing times",
      classifiers=["Programming Language :: Python :: 3"],
      )
