"""File containing functions for reconstruction of the pair distribution
function."""

import os.path
import numpy as np

from .fitfuncs import FOURIERFITFUNCS, PARAMETERFITFUNCS

__all__ = ["loadParameterFile", "reconstruct_mgU2prime", "getU2prime"]

# -- Internal constants for loading the default parameter file --

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
FITPARAMSFILE = "fit_parameters.csv"
DEFAULTPATH = os.path.join(SCRIPTDIR, FITPARAMSFILE)
DEFAULTPARAMETERS = None # Lazy loading: will be set when needed

# -- Simulation constants --

EPSILON, SIGMA = 1,1

# -- File input --

def loadParameterFile(filepath):
    r"""Load a parameter CSV file from a given path.

    Parameters
    ----------
    path: str
        Path to CSV file.

    Returns
    -------
    params_f11, params_f22: dict
        Dictionaries containing all fit parameters for both f_11 and f_22
        Fourier coefficients.
    """
    params_f11 = {}
    params_f22 = {}
    # Open file with fit parameters and read line by line every q_i,j
    # or v,w resp.
    with open(filepath) as f:
        for line in f:
            # Extract cells
            cells = line.split(",")
            # Skip lines without label
            if cells[0] == "":
                continue
            # Skip invalid labels
            if len(cells[0]) != 7:
                print("Invalid Fourier coefficient label: {}".format(cells[0]))
                continue
            coefftype = cells[0][1:4] # Either "^11" or "^22"
            k,l = map(int, cells[0][5:])
            # Select correct dict
            params = params_f11 if coefftype == "^11" else params_f22
            # Filter empty cells
            cells = list(filter(lambda x: x != "", cells))
            if not (k,l) in params:
                params[(k,l)] = []
            params[(k,l)].append(list(map(float, cells[2:]))) # Extract data
        return params_f11, params_f22

# Internally used
def loadDefaultParameterFile():
    global DEFAULTPARAMETERS
    if not DEFAULTPARAMETERS:
        DEFAULTPARAMETERS = loadParameterFile(DEFAULTPATH)
    return DEFAULTPARAMETERS

# -- Reconstruction functions --

def reconstruct_mgU2prime(
        r, phi1, phi2, phi0, Pe, params_f11=None, params_f22=None):
    r"""Returns an approximation for -gU'_2 in a given range of particle
    distances and positional and orientational angles.

    Parameters
    ----------
    r: float or array_like
        Distance(s) at which -gU'_2 will be calculated
    phi1, phi2: float or array_like
        Positional and orientational angles at which -gU'_2 will be calculated
    phi0, Pe: float
        Packing density and Peclet number for which -gU'_2 will be calculated
    params_f11, params_f22: dict
        Parameter dictionary containing all fit parameters necessary for
        reconstruction. If not set, the included default values will be used.
    """
    # If params are not set: load default
    if not params_f11 or not params_f22:
        params = loadDefaultParameterFile()
        if not params_f11:
            params_f11 = params[0]
        if not params_f22:
            params_f22 = params[1]
    # r must have dimension 1
    if len(np.shape(r)) == 0:
        r = np.array([r])
    # Allocate array for -gU'_2
    gU2 = np.zeros((r.shape[0],
        1 if len(np.shape(phi1)) == 0 else phi1.shape[0],
        1 if len(np.shape(phi1)) == 0 else phi1.shape[1]))
    # Calculate all contributions for k,l between 0 and 2
    for k in range(3):
        for l in range(3):
            # Start with cos*cos for f^11 coefficient
            fourierfunc = np.cos
            for params in (
                    [params_f11] if k == 0 or l == 0 else
                    [params_f11, params_f22]):
                # Calculate fit parameters f_0, mu, sigma, lambda and r_1, r_2
                # (if necessary)
                fitparams = list(map(lambda x : x[0]((phi0, Pe), *x[1]),
                    zip(PARAMETERFITFUNCS[(k,l)], params[(k,l)])))
                # Calculate Fourier coefficient
                fouriercoeff = FOURIERFITFUNCS[(k,l)](r, *fitparams)
                # Add contribution to -gU'_2
                gU2 += fouriercoeff[:, None, None] * \
                fourierfunc(phi1*k) * fourierfunc(phi2*l)
                # Take sin*sin next time (for f^22 coefficient)
                fourierfunc = np.sin
    return gU2

def getU2prime(r):
    r"""Calculate derivative of potential U_2 with respect to r.

    The potential depends on the particle diameter sigma and the Lennard-Jones
    energy epsilon. This function uses values consistent with the simulations
    described in the accompanying article.

    Parameters
    ----------
    r: float or array_like
        Distance(s) at which the derivative of the interaction potential will be
        calculated.
    """
    return 24*EPSILON/SIGMA * (-2*(SIGMA/r)**13 + (SIGMA/r)**7)
