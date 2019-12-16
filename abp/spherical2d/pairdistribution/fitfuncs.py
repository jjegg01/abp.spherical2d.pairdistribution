"""Functions for the fitting procedures of the Fourier coefficients of -gU'_2.

For the functions of the first fitting procedure the order of arguments matches
the one shown in the row headers for each Fourier coefficient in fit_parameters.csv.
Likewise for the functions of the second fitting procedure the order of
arguments matches the order of arguments shown in the header line in
fit_parameters.csv."""

from numpy import sqrt, exp
from scipy.special import erfc #pylint: disable=no-name-in-module

# -- Functions for first fitting procedure --

def EMG(r, mu, sig, la):
    """Exponentially modified Gaussian distribution."""
    return la/2 * exp(la/2 * (la*sig**2 - 2*(r-mu))) * \
        erfc((la*sig**2 - (r-mu))/sqrt(2)/sig)

def FourierFit0(r,a,mu,sig,la):
    """Function for first fitting procedure with one fixed root."""
    return a*EMG(r,mu,sig,la) * (2**(1/6) - r)

def FourierFit1(r,a,mu,sig,la,b):
    """Function for first fitting procedure with one fixed and one variable
    root."""
    return FourierFit0(r,a,mu,sig,la) * (r-b)

def FourierFit2(r,a,mu,sig,la,b,c):
    """Function for first fitting procedure with one fixed and two variable
    roots."""
    return FourierFit1(r,a,mu,sig,la,b) * (r-c)

FOURIERFITFUNCS = {
    (0,0) : FourierFit0,
    (1,0) : FourierFit1,
    (2,0) : FourierFit2,
    (0,1) : FourierFit1,
    (1,1) : FourierFit1,
    (2,1) : FourierFit2,
    (0,2) : FourierFit1,
    (1,2) : FourierFit2,
    (2,2) : FourierFit1
}
"""Dictionary mapping the index tuple (k,l) to the corresponding fit function"""

# -- Functions for second fitting procedure --

def h2(x, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15):
    """Function h_2 for second fitting procedure."""
    dens, Pe = x
    return q1/Pe + q2/sqrt(Pe) + q3 + q4*sqrt(Pe) + q5*Pe + \
          (q6/Pe + q7/sqrt(Pe) + q8 + q9*sqrt(Pe) + q10*Pe)*dens + \
          (q11/Pe + q12/sqrt(Pe) + q13 + q14*sqrt(Pe) + q15*Pe)*dens**2

def h2_0(x, *args):
    """Function h_{2,0} for second fitting procedure."""
    #pylint: disable=no-value-for-parameter,unbalanced-tuple-unpacking
    if len(args) != 17:
        raise TypeError("h2_0 takes exactly 17 arguments ({} given)".format(len(args)))
    dens, _ = x
    v, w = args[-2:]
    return h2(x, *(args[:-2])) + v * exp(w*dens)

def h3(x, *args):
    """Function h_3 for second fitting procedure."""
    #pylint: disable=no-value-for-parameter,unbalanced-tuple-unpacking
    if len(args) != 20:
        raise TypeError("h3 takes exactly 20 arguments ({} given)".format(len(args)))
    dens, Pe = x
    q16, q17, q18, q19, q20 = args[-5:]
    return h2(x, *(args[:-5])) + \
        (q16/Pe + q17/sqrt(Pe) + q18 + q19*sqrt(Pe) + q20*Pe)*dens**3

def h2_1(x, *args):
    """Function h_{2,1} for second fitting procedure."""
    #pylint: disable=no-value-for-parameter,unbalanced-tuple-unpacking
    if len(args) != 17:
        raise TypeError("h2_1 takes exactly 17 arguments ({} given)".format(len(args)))
    dens, Pe = x
    v, w = args[-2:]
    return h2(x, *(args[:-2])) + v * exp(w*dens)/Pe

#
PARAMETERFITFUNCS = {
    (0,0) : [h2_1, h2, h2, h2_0],
    (0,1) : [h3, h2, h2, h2_0, h3],
    (0,2) : [h3, h2, h2, h2_0, h3],
    (1,0) : [h3, h2, h2, h2_0, h3],
    (1,1) : [h3, h2, h2, h2_0, h3],
    (1,2) : [h3, h2, h2, h2_0, h3, h3],
    (2,0) : [h3, h2, h2, h2_0, h3, h3],
    (2,1) : [h3, h2, h2, h2_0, h3, h3],
    (2,2) : [h3, h2, h2, h2_0, h3],
}
"""Dictionary mapping the index tuple (k,l) to an array of fit functions for
each fit parameter of FourierFit{0,1,2}."""
