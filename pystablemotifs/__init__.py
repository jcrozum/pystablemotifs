try:
    import pyboolnet
except ImportError as e:
    print(e)
    print(
        """
pyboolnet can be installed via
pip install git+https://github.com/hklarner/pyboolnet
or, for a specific version, via, e.g.,
pip install git+https://github.com/hklarner/pyboolnet@3.0.9
"""
    )
    raise

from pystablemotifs import (
    drivers,
    export,
    format,
    random_boolean_networks,
    reduction,
    restrict_space,
    succession,
    time_reversal,
)
from pystablemotifs.Attractor import Attractor
from pystablemotifs.AttractorRepertoire import AttractorRepertoire

__version__ = "3.0.6"
