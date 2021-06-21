
try:
    import PyBoolNet
except ImportError as e:
    print(e)
    print("""
PyBoolNet can be installed via
pip install git+https://github.com/hklarner/PyBoolNet
or, for a specific version, via
pip install git+https://github.com/hklarner/PyBoolNet@2.30.0
""")
    raise

from PyStableMotifs import Format
from PyStableMotifs import DomainOfInfluence
from PyStableMotifs import Reduction
from PyStableMotifs import TimeReversal
from PyStableMotifs import RestrictSpace
from PyStableMotifs import Succession
from PyStableMotifs import RandomBooleanNetworks
from PyStableMotifs.Attractor import Attractor
from PyStableMotifs.AttractorRepertoire import AttractorRepertoire
from PyStableMotifs import Export

__version__ = '2.1.0'
