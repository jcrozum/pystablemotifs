
try:
    import pyboolnet
except ImportError as e:
    print(e)
    print("""
pyboolnet can be installed via
pip install git+https://github.com/hklarner/pyboolnet
or, for a specific version, via
pip install git+https://github.com/hklarner/pyboolnet@3.0.5
""")
    raise

from pystablemotifs import format
from pystablemotifs import drivers
from pystablemotifs import reduction
from pystablemotifs import time_reversal
from pystablemotifs import restrict_space
from pystablemotifs import succession
from pystablemotifs import random_boolean_networks
from pystablemotifs.Attractor import Attractor
from pystablemotifs.AttractorRepertoire import AttractorRepertoire
from pystablemotifs import export

__version__ = '2.1.0'
