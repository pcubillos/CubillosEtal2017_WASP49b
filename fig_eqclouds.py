#! /usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d as gaussf

def getprofiles(file):
  """
  Read Helios temperature-pressure profiles.
  """
  with open(file) as f:
    lines = f.readlines()
    nlines  = len(lines)
    nlayers = int(lines[-1].split()[0]) + 1

    temp  = np.zeros(nlayers, np.double) # Kelvin
    press = np.zeros(nlayers, np.double) # Barye

    for i in np.arange(nlayers):
      temp[i], press[i] = lines[i+nlines-nlayers].split()[1:]

  return temp, press

# HELIOS TP files:
helios = ["./inputs/helios/tp_w49b_100x_solar.dat",
          "./inputs/helios/tp_cloud_higher_stronger.dat",
          "./inputs/helios/tp_w49b_0.1x_solar.dat"]

# Equilibrium condensation-curve files:
eqclouds = ["./inputs/eqclouds/Twasp49blarge.dat",
            "./inputs/eqclouds/Twasp49bsolar.dat",
            "./inputs/eqclouds/Twasp49bsmall.dat"] 

condensates = np.array([
    r"${\rm Cr}$",      r"${\rm Fe}$",   r"${\rm Mg_2SiO_4}$",
    r"${\rm MgSiO_3}$", r"${\rm Mn}$",   r"${\rm MnS}$",
    r"${\rm Na_2S}$",   r"${\rm ZnS}$",  r"${\rm KCl}$",
    r"${\rm NH_4SH}$",  r"${\rm H_2S}$", r"${\rm NH_4H_2PO_4}$"])

with open(eqclouds[0]) as f:
  lines = f.readlines()

nlayers = len(lines) - 13
ncond   = len(condensates)
nz      = len(eqclouds)

press = np.zeros((nz, nlayers),        np.double)
cond  = np.zeros((nz, ncond, nlayers), np.double)

htemp = []
for k in np.arange(nz):
  # Read HELIOS file:
  temp, hpress = getprofiles(helios[k])
  htemp.append(temp)
  hpress /= 1e6  # in bar
  # Read eq. cloud file:
  with open(eqclouds[k]) as f:
    lines = f.readlines()
  for i in np.arange(nlayers):
    info = lines[i+13].split()
    press[k,i]  = info[0]
    cond[k,:,i] = info[1:]

# Connect to thermosphere:
a, b = -1100.0, -17600
thermo = a*np.log(hpress) + b  # Temperature slope as estimade from HD runs
for k in np.arange(nz):
  ithermo = thermo > htemp[k]
  htemp[k][ithermo] = thermo[ithermo]
  # Round the intersection:
  iround = np.where(ithermo)[0][0] - 3
  htemp[k][iround:] = gaussf(htemp[k][iround:], 2)

# Read transit pressure HPD boundaries:
with open("run02/pt_boundaries.dat", "r") as f:
  info = f.readlines()

HPD = np.zeros((3,2,2), np.double)
for i in np.arange(nz):
  hpd = info[i+1].split()
  HPD[i,0] = hpd[0:2]
  HPD[i,1] = hpd[2:4]
HPD = 10**HPD

# Plot setup:
yran = 1e2, 0.3e-8
xran = 0, 2990
hw = 100
matplotlib.rcParams.update({'xtick.minor.size':0.0})
matplotlib.rcParams.update({'xtick.labelsize':10})
# colors
c = np.array(['aqua',  'orange', 'brown',     'limegreen', 'k',
              'gold',  'peru',   'darkgreen', 'r',
              'olive', 'orchid', 'b'])

csort = [ 2, 3,  1, 0, 5,
          4, 6,  7, 8,
         11, 9, 10]
condensates = condensates[csort]

# Line style:
d = [(), (), (), (), (),                # Solid
     (5,2), (5,2), (5,2), (5,2),        # Short dash
     (5,2,2,2), (5,2,2,2), (5,2,2,2), ] # long dash

# Metallicity label:
xtext, ytext = 2500, 9e-9
text = [r"$100\times\,{\rm solar}$",
        r"$1.0\times\,{\rm solar}$",
        r"$0.1\times\,{\rm solar}$"]

# The plot:
plt.figure(-5, (8.5, 5))
plt.clf()
plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.96,
                    hspace=0.13, wspace=0.01)
k=0
ax=plt.subplot(1, 3, k+1)
plt.axhspan(HPD[k][1][0], HPD[k][1][1], color="0.925", zorder=0)  # 95%
plt.axhspan(HPD[k][0][0], HPD[k][0][1], color="0.8",   zorder=0)  # 68%
plt.fill_betweenx(hpress, htemp[k]-hw, htemp[k]+hw,
                  color="cornflowerblue", zorder=1)
for i in np.arange(ncond):
  plt.semilogy(cond[k,csort[i]], press[k], lw=1.0, color=c[i], dashes=d[i])
plt.ylabel(r"${\rm Pressure\ \ (bar)}$")
plt.xlabel(r"${\rm Temperature\ \ (K)}$")
plt.text(xtext, ytext, text[k], fontsize=12, ha="right")
plt.ylim(yran)
plt.xlim(xran)

k += 1
ax=plt.subplot(1, 3, k+1)
plt.axhspan(HPD[k][1][0], HPD[k][1][1], color="0.925", zorder=0)
plt.axhspan(HPD[k][0][0], HPD[k][0][1], color="0.8",   zorder=0)
plt.fill_betweenx(hpress, htemp[k]-hw, htemp[k]+hw,
                  color="cornflowerblue", zorder=1)
for i in np.arange(ncond):
  plt.semilogy(cond[k,csort[i]], press[k], lw=1.0, color=c[i], dashes=d[i])
plt.xlabel(r"${\rm Temperature\ \ (K)}$")
plt.text(xtext, ytext, text[k], fontsize=12, ha="right")
plt.ylim(yran)
plt.xlim(xran)
ax.set_yticklabels([])

k += 1
ax=plt.subplot(1, 3, k+1)
plt.axhspan(HPD[k][1][0], HPD[k][1][1], color="0.925", zorder=0)
plt.axhspan(HPD[k][0][0], HPD[k][0][1], color="0.8",   zorder=0)
plt.fill_betweenx(hpress, htemp[k]-hw, htemp[k]+hw,
                  color="cornflowerblue", zorder=1)
for i in np.arange(ncond):
  plt.semilogy(cond[k,csort[i]], press[k], lw=1.0, color=c[i],
               label=condensates[i], dashes=d[i])
plt.xlabel(r"${\rm Temperature\ \ (K)}$")
plt.text(xtext, ytext, text[k], fontsize=12, ha="right")
plt.ylim(yran)
plt.xlim(xran)
plt.legend(loc="upper right", fontsize=7.5, bbox_to_anchor=(0.99, 0.57),
           bbox_transform=plt.gcf().transFigure)
ax.set_yticklabels([])
plt.savefig("./plots/equilibrium_clouds.ps")
