import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import abp.spherical2d.pairdistribution as abp2dpd

# -- Parse args --

parser = argparse.ArgumentParser(
    description=r"""Display the pair-distribution function for a given particle
        distance, Peclet number and packing density. The default values are the
        same as in Fig. 3 of the accompanying article by J. Jeggle,
        J. Stenhammar and R. Wittkowski.""")
parser.add_argument(
    "-r", metavar="dist", dest="dist", type=float, default=1.0,
    help="Particle distance in multiples of sigma (default: 1)")
parser.add_argument(
    "-d", metavar="phi0", dest="phi0", type=float, default=0.2,
    help="Packing density (default: 0.2)")
parser.add_argument(
    "-p", metavar="peclet", dest="peclet", type=float, default=50.0,
    help="Peclet number (default: 50)")

args = parser.parse_args()

# Validate args

r_min = 0.7775
r_max = 2**(1/6)
if args.dist < r_min or args.dist > r_max:
    print("Warning: Distance is outside of approximation bounds!")

if args.peclet < 0:
    print("Warning: Unphysical argument for Peclet number")
if args.phi0 < 0 or args.phi0 > 1:
    print("Warning: Unphysical argument for packing density")

# -- Calculate pair-distribution function --

# Generate arrays for r, phi1 and phi2
phi1 = np.linspace(-np.pi,np.pi,180,endpoint=False)
phi2 = np.linspace(0,2*np.pi,180,endpoint=False)
phi1,phi2 = np.meshgrid(phi1, phi2)
r = args.dist # Just take a single distance

# Calculate -gU'_2
gU2 = abp2dpd.reconstruct_mgU2prime(r, phi1, phi2, args.phi0, args.peclet)[0]

# Divide by U'_2 to obtain the pair-distribution function g
g = -gU2/abp2dpd.getU2prime(args.dist)

# -- Plotting code --

fig, ax = plt.subplots(1)

test = np.zeros((11,100))
cax = ax.imshow(g.T, cmap="inferno", origin="lower",
extent=(0,g.shape[0],0,g.shape[0]))
cbar = fig.colorbar(cax)

cbar.set_label("$g$")

ax.set_xlabel(r"$\phi_1$")
ax.set_ylabel(r"$\phi_2$")

ax.set_xticks([0,g.shape[0]//2,g.shape[0]])
ax.set_xticklabels([r"$-\pi$",r"0",r"$\pi$"])
ax.xaxis.set_minor_locator(MultipleLocator(g.shape[0]//4))

ax.set_yticks([0,g.shape[0]//2,g.shape[0]])
ax.set_yticklabels(["0",r"$\pi$",r"$2\pi$"])
ax.yaxis.set_minor_locator(MultipleLocator(g.shape[0]//4))

plt.show()
