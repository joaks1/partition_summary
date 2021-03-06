from setuptools import setup, find_packages
import sys, os

SCRIPT_SUBPATHS = [
    ['partitionpy', 'partition_sum.py'],
 ]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]

version = '0.1'

setup(name='partitionpy',
      version=version,
      description="Python package for manipulating and summarizing posteriors of DNA alignment partitions produced by Bayesian Dirichlet process models of across-site rate heterogeneity",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='phylogenetics DPP evolution phylogeny genetics',
      author='Jamie R. Oaks and Mark T. Holder',
      author_email='joaks1@gmail.com and mtholder@ku.edu',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts = SCRIPTS,
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          'munkres'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
