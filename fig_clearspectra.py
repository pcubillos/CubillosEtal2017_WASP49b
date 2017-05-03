#! /usr/bin/env python

import sys, os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
import scipy.integrate   as si
import scipy.interpolate as sip
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ion()

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc


# Plot clear spectra for each metallicity:
pyrat = pb.pbay.run("clear_spectra.cfg")

# Read atmospheric models:
mol, press, temp, abun1  = pb.atmosphere.readatm("WASP-49b_000.1xsolar.atm")
mol, press, temp, abun2  = pb.atmosphere.readatm("WASP-49b_001.0xsolar.atm")
mol, press, temp, abun3  = pb.atmosphere.readatm("WASP-49b_100.0xsolar.atm")

# Compute transmission spectra for each atmosphere:
pyrat = pb.pyrat.run(pyrat, [temp, abun1])
rp1  = np.sqrt(pyrat.spec.spectrum)*pyrat.phy.rstar  # Transit radius
rad1 = pyrat.atm.radius                  # Atmospheric radius profile

pyrat = pb.pyrat.run(pyrat, [temp, abun2])
rp2  = np.sqrt(pyrat.spec.spectrum)*pyrat.phy.rstar
rad2 = pyrat.atm.radius

pyrat = pb.pyrat.run(pyrat, [temp, abun3])
rp3  = np.sqrt(pyrat.spec.spectrum)*pyrat.phy.rstar
rad3 = pyrat.atm.radius


# Pressure at transit radius as function of wavelength:
sigma = 3.0  # Gauss-convolve for better-looking plots
p   = sip.interp1d(rad1[::-1], press[::-1])
pt1 = p(gaussf(rp1, sigma))
p   = sip.interp1d(rad2[::-1], press[::-1])
pt2 = p(gaussf(rp2, sigma))
p   = sip.interp1d(rad3[::-1], press[::-1])
pt3 = p(gaussf(rp3, sigma))

# Photospheric pressure modulation spectrum:
lw = 1.5
plt.figure(-25)
plt.clf()
ax = plt.subplot(111)
plt.semilogy(1e4/pyrat.spec.wn, pt3, lw=lw, color="orange",
             label=r"$100\times\,{\rm solar}$")
plt.semilogy(1e4/pyrat.spec.wn, pt2, lw=lw, color="sienna",
             label=r"$1.0\times\,{\rm solar}$")
plt.semilogy(1e4/pyrat.spec.wn, pt1, lw=lw, color="k",
             label=r"$0.1\times\,{\rm solar}$")
plt.axvspan(0.74, 1.01, color="0.9")

plt.xlim(0.5, 1.2)
plt.ylim(3, 3e-7)
plt.legend(loc="upper right",  fontsize=15)
plt.ylabel(r"$\rm Pressure\ \ (bar)$",  fontsize=16)
plt.xlabel(r"$\rm Wavelength\ \ (um)$", fontsize=16)
plt.savefig("../plots/WASP49b_clear_spectra.ps")

