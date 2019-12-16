"""Module for approximating the negative product of the pair-distribution function and
the derivative of the interaction potential for a two-dimensional suspension
of spherical active Brownian particles."""

__author__ = "Julian Jeggle, Raphael Wittkowski"
__copyright__ = "Copyright (C) 2019 Julian Jeggle, Raphael Wittkowski"
__license__ = "MIT"
__version__ = "1.0"

from ._reconstruct import *

# Trick pdoc3 into documenting submodule contents as members of this module
from ._reconstruct import __all__ as reconstruct_all
__all__ = reconstruct_all
