#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gaussf

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.utils as mu

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc


# Load HARPS data from Wittenbach et al. (2017):
with open("../inputs/WASP49b_HARPS_wyttenbach.rdb") as f:
  lines = f.readlines()

header = 2
nwave = len(lines) - header

wlength = np.zeros(nwave, np.double)
rtilde  = np.zeros(nwave, np.double)
uncert  = np.zeros(nwave, np.double)

for i in np.arange(nwave):
  wlength[i], rtilde[i], uncert[i] = lines[i+header].split()

binsize = 15
binr, binu, binwl = mu.binarray(rtilde, uncert, wlength, binsize)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# High resolution RT model:
pyrat = pb.pyrat.init("spectra_hires.cfg")

q  = pyrat.atm.q
wl = 1e4/pyrat.spec.wn
# Doppler-shifted wl:
vel = -41.726 * pc.km
dwl = wl * (1-vel/pc.c) / np.sqrt(1-(vel/pc.c)**2)

# Mean wavelength sampling rate in microns:
deltawl = np.abs(np.mean(np.ediff1d(wl)))
wl0 = np.mean(wl)
# Set resolution R=115000 at lambda=wl0
R = 115000.0
dlambda = wl0/R  # lambda/R
# Standard-deviation corresponding to instrumental resolution:
sigma = (dlambda/deltawl)/2.355

mol, press, temp, q1 = pb.atmosphere.readatm("WASP-49b_000.1xsolar_hires.atm")
mol, press, temp, q2 = pb.atmosphere.readatm("WASP-49b_001.0xsolar_hires.atm")
mol, press, temp, q3 = pb.atmosphere.readatm("WASP-49b_100.0xsolar_hires.atm")


# Radius-fgray values from pair correlation plots:
pyrat.haze.model[0].pars[0] = 1.21
pyrat.phy.rplanet = 1.164*pc.rjup
pyrat = pb.pyrat.run(pyrat, [temp, q1])
spectrum1 = pyrat.spec.spectrum
M1 = (1-spectrum1) / (1-spectrum1[0])

pyrat.haze.model[0].pars[0] = 1.6
pyrat.phy.rplanet = 1.16*pc.rjup
pyrat = pb.pyrat.run(pyrat, [temp, q2])
spectrum2 = pyrat.spec.spectrum
M2 = (1-spectrum2) / (1-spectrum2[0])

pyrat.haze.model[0].pars[0] = 3.63
pyrat.phy.rplanet = 1.156*pc.rjup
pyrat = pb.pyrat.run(pyrat, [temp, q3])
spectrum3 = pyrat.spec.spectrum
M3 = (1-spectrum3) / (1-spectrum3[0])


# The plot:
plt.figure(5, (8.5,5))
plt.clf()
ax = plt.axes([0.11, 0.2, 0.6, 0.55])
# The models:
plt.plot(1e4*dwl, gaussf(M3, sigma), "-", lw=1.25, zorder=0,
         label=r"$100\times\,{\rm solar}$", color='orange')
plt.plot(1e4*dwl, gaussf(M2, sigma), "-", lw=1.25, zorder=0,
         label=r"$1.0\times\,{\rm solar}$", color='sienna')
plt.plot(1e4*dwl, gaussf(M1, sigma), "-", lw=1.25, zorder=0,
         label=r"$0.1\times\ {\rm solar}$", color='k')
# The data:
plt.errorbar(binwl, binr, binu, fmt="o", color="0.6", ms=4, zorder=1,
             label=r"$\rm Wyttenbach\ et\ al.\ (2017)$")
ax.set_xticks([5890, 5892, 5894, 5896, 5898, 5900])
ax.set_xticklabels(['5890', '5892', '5894', '5896', '5898', '5900'])
plt.xlim(5889, 5899)
plt.ylim(0.972, 1.0095)
plt.ylabel(r"$\tilde{\mathfrak{R}}$", fontsize=14)
plt.xlabel(r"$\rm Wavelength\ (um)$", fontsize=14)
plt.legend(loc=(0.25, 0.04), fontsize=11)
# Temperature profile:
ax = plt.axes([0.815, 0.2, 0.15, 0.55])
plt.semilogy(temp, pyrat.atm.press/pc.bar, "r-", lw=1.5)
plt.xlim(900, 3100)
plt.ylim(np.amax(pyrat.atm.press/pc.bar), np.amin(pyrat.atm.press/pc.bar))
plt.xticks([1000, 2000, 3000])
plt.xlabel(r"$\rm Temperature\ (K)$", fontsize=14)
plt.ylabel(r"$\rm Pressure\ (bar)$", fontsize=14)
plt.savefig("../plots/WASP49b_hires.ps")

