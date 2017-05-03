#! /usr/bin/env python

import sys
import matplotlib
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

sys.path.append("./pyratbay")
import pyratbay.constants as pc


def getprofiles(file):
  with open(file) as f:
    lines = f.readlines()
    nlines = len(lines)
    nlayers = int(lines[-1].split()[0]) + 1
    temp  = np.zeros(nlayers, np.double) # Kelvin
    press = np.zeros(nlayers, np.double) # Barye
    for i in np.arange(nlayers):
      temp[i], press[i] = lines[i+nlines-nlayers].split()[1:]
  return temp, press


# Read HELIOS radiative-equilibrium profiles:
temp = []
t, press = getprofiles("./inputs/helios/tp_clearsky.dat")
temp.append(t)
t, press = getprofiles("./inputs/helios/tp_cloud_higher_stronger.dat")
temp.append(t)
t, press = getprofiles("./inputs/helios/tp_cloud_higher_weaker.dat")
temp.append(t)
t, press = getprofiles("./inputs/helios/tp_cloud_lower_stronger.dat")
temp.append(t)
t, press = getprofiles("./inputs/helios/tp_cloud_lower_weaker.dat")
temp.append(t)

# Pressure in bar:
press /= 1e6

matplotlib.rcParams.update({'ytick.minor.size':0.0})
fs = 17
lw = 2

# The plot:
plt.figure(-11, (7,6))
plt.clf()
ax = plt.subplot(111)
plt.semilogy(temp[0], press, lw=lw, color="red",       dashes=(8,2,3,2),
             label=r"$\rm Clear$")
plt.semilogy(temp[2], press, lw=lw, color="navy",
             label=r"$\rm Cloud\ at\ 10^{-6}\ bar,\ weaker$")
plt.semilogy(temp[1], press, lw=lw, color="navy",      dashes=(7,3),
             label=r"$\rm Cloud\ at\ 10^{-6}\ bar,\ stronger$")
plt.semilogy(temp[4], press, lw=lw, color="limegreen",
             label=r"$\rm Cloud\ at\ 10^{-5}\ bar,\ weaker$")
plt.semilogy(temp[3], press, lw=lw, color="limegreen", dashes=(7,3),
             label=r"$\rm Cloud\ at\ 10^{-5}\ bar,\ stronger$")
plt.ylim(np.amax(press), np.amin(press))
plt.ylim(1e2, 1e-8)
plt.ylabel(r"${\rm Pressure\ \ (bar)}$",  fontsize=fs)
plt.xlabel(r"${\rm Temperature\ \ (K)}$", fontsize=fs)
plt.legend(loc="upper right")
plt.savefig("./plots/radeq_temperature.ps")

