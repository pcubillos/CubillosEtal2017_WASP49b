#! /usr/bin/env python

import os, sys
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

sys.path.append("./pyratbay")
import pyratbay.constants as pc


# Read hydrodynamic profiles:
with open("./inputs/hydro/Pressure.dat") as f:
  lines = f.readlines()
with open("./inputs/hydro/Temperature.dat") as f:
  tlines = f.readlines()
with open("./inputs/hydro/H_density.dat") as f:
  Hlines = f.readlines()
with open("./inputs/hydro/H2_density.dat") as f:
  H2lines = f.readlines()

nlayers = len(lines)
temp   = np.zeros(nlayers, np.double)
H2dens = np.zeros(nlayers, np.double)
Hdens  = np.zeros(nlayers, np.double)
rad    = np.zeros(nlayers, np.double)
press  = np.zeros(nlayers, np.double)
for i in np.arange(nlayers):
  rad[i], press[i] = lines[i].split()
  temp  [i] = tlines[i].split()[1] 
  Hdens [i] = Hlines[i].split()[1] 
  H2dens[i] = H2lines[i].split()[1] 

n  = 1e-1 * press / (sc.k * temp)

# The plot:
plt.figure(-1, (7,6))
plt.clf()
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95,
                    wspace=0.05)
ax = plt.subplot(131)
plt.semilogy(rad/pc.rjup,  press,  "r-", lw=2,
             label=r"$1\!\times\,F_{\rm XUV}$", zorder=4)
plt.ylabel(r"${\rm Pressure\ \ (bar)}$",       fontsize=14)
plt.xlabel(r"${\rm Radius}\ \ (R_{\rm jup})$", fontsize=14)
plt.ylim(np.amax(press), np.amin(press))
plt.xticks([1, 2, 3, 4.0])
ax.set_xticklabels(["1.0", "2.0", "3.0", "4.0"])

plt.subplot(132)
plt.semilogy(temp,  press,  "r-", lw=2)
plt.xlabel(r"${\rm Temperature\ \ (K)}$", fontsize=14)
plt.ylim(np.amax(press), np.amin(press))
plt.xticks([1000, 2000, 3000, 4000])
plt.xlim(500, 4500)
a = plt.yticks(visible=False)

plt.subplot(133)
plt.loglog(Hdens /n, press,"-",lw=2,color="limegreen",  label=r"$\rm H$")
plt.loglog(H2dens/n, press,"-",lw=2,color="navy", label=r"${\rm H}_2$")
plt.xlabel(r"${\rm Mixing\ ratio}$", fontsize=14)
plt.ylim(np.amax(press), np.amin(press))
plt.xticks([1e-9, 1e-6, 1e-3, 1])
a = plt.yticks(visible=False)
plt.xlim(1e-7, 10.0)
plt.legend(loc='best')
plt.savefig("./plots/WASP49b_hydro_profiles.ps")

