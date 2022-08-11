import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.0.1'
PACKAGE_NAME = 'ReliabilityDiagram'
AUTHOR = 'Naveen Goutham, Camille Le Coz'
AUTHOR_EMAIL = 'naveen.goutham@outlook.com'
URL = 'https://github.com/clecoz/ReliabilityDiagram'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'A package to compute contingency table, forecast attributes, and reliability diagram for forecasts.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = ['numpy','statsmodels','matplotlib']

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      keywords=['reliability table','contingency table','reliability diagram','attributes diagram','reliability','resolution','forecast attributes','ensemble forecast','forecast','climatology','python','matplotlib'],
      packages=find_packages()
      )
