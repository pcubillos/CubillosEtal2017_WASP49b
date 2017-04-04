#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

sys.path.append("../")
import hists as h

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed as mc3
sys.path.append("../pyratbay")
import pyratbay as pb
pc  = pb.constants
atm = pb.atmosphere

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def readmolfile(molfile):
  """
  Document me!
  """
  with open(molfile, "r") as molfile:
    # Skip comment and blank lines:
    line = molfile.readline().strip()
    while line == '' or line.startswith('#'):
      line = molfile.readline().strip()

    molID  = [] # Molecule ID
    symbol = [] # Molecule symbol
    mass   = [] # Molecule mass
    diam   = [] # Molecule diameter
    # Read Molecular values:
    while line != '' and not line.startswith('#'):  # Start reading species
      molinfo = line.split()
      # Extract info:
      molID .append(  int(molinfo[0]))
      symbol.append(      molinfo[1] )
      mass  .append(float(molinfo[2]))
      diam  .append(float(molinfo[3]))
      line = molfile.readline().strip()  # Read next line

  molID  = np.asarray(molID)
  symbol = np.asarray(symbol)
  mass   = np.asarray(mass)
  diam   = np.asarray(diam)

  return molID, symbol, mass, diam

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ID, symbol, mass, diam = readmolfile("../pyratbay/inputs/molecules.dat")

p0       = 0.1   * pc.bar   # Reference pressure level
rtransit = 1.198 * pc.rjup  # Transit radius (white light)
mplanet  = 0.396 * pc.mjup  # Planetary mass
rstar    = 1.038 * pc.rsun  # Stellar radius

atmfile = ["../run01/WASP-49b_100.0xsolar.atm",
           "../run01/WASP-49b_001.0xsolar.atm",
           "../run01/WASP-49b_000.1xsolar.atm"]

postfiles = ["./MCMC_WASP-49b_000.1xsolar_Tprior.npz",
             "./MCMC_WASP-49b_000.1xsolar_Tprior.npz",
             "./MCMC_WASP-49b_000.1xsolar_Tprior.npz"]

name = ["100.0xsolar", "001.0xsolar", "000.1xsolar"]

parname = [r"$\log_{10}(p_{\rm T})\ \ ({\rm bar})$",
           r"$T\ \ ({\rm K})$",
           r"$R_{0.1\,{\rm bar}}\ \ (R_{\rm jup})$",
           r"$\log_{10}(f_{\rm gray})$"]

burn = 10000
Zp, Zr = [], []
for j in np.arange(3):
  # Read MCMC results:
  post = np.load(postfiles[j])["Z"]
  temp = post[burn:,0]
  rad  = post[burn:,1] * pc.km
  cs   = post[burn:,2]
  
  # Read atmosphere:
  species, press, t, q = atm.readatm(atmfile[j])
  press *= pc.bar
  logp = np.log10(press)
  # Mean molecular mass
  mspec = []
  for spec in species:
    mspec.append(mass[symbol==spec][0])
  mspec = np.asarray(mspec)
  mu      = np.sum(q*mspec, axis=1) # Mean molecular mass profile
  
  # Unique posterior values:
  utemp, uind, uinv = np.unique(temp, return_index=True, return_inverse=True)
  urad   = rad[uind]
  upress = np.zeros(len(utemp), np.double)
  
  nsmall = 0
  for i in np.arange(len(utemp)):
    # Compute radius atmospheric profile:
    radius = atm.hydro_m(press, utemp[i], mu, mplanet, p0, urad[i])
    if radius[0] < rtransit:
      nsmall += 1
      continue
    # Interpolate pressure-radius:
    f = si.interp1d(radius, logp)
    # Get pressure of the photospheric radius (transit radius):
    upress[i] = 10**f(rtransit) / pc.bar
    if i % (len(utemp)/8) == 0:
      print("{:6.2f}%  ({:d})".format(i*100.0/len(utemp), i))
  # Posterior distribution of the photospheric pressure of the planet:
  pressure = upress[uinv]
  
  # Remove instances where the planet atmosphere does not reach rtransit:
  good = pressure > 0
  pressure    = pressure[good]
  temperature = temp[good]
  radius      = rad[good] / pc.rjup
  graycs      = cs[good]
 
  Z = np.vstack((np.log10(pressure), temperature, radius, graycs)).T
  Zp.append(np.log10(pressure))
  Zr.append(radius)
  # Figure 4:
  mc3.plots.pairwise(Z, parname=parname, fignum=j+10,
           savefile="../plots/pair_{}.ps".format(name[j]))
  mc3.plots.histogram(Z, parname=parname, percentile=0.683, fignum=20+j)

xran = [[-8, -1.0], [1.06, 1.205]]
pname = [r"$\log_{10}(p_{\rm T})\ \ ({\rm bar})$",
         r"$R_{0.1\,{\rm bar}}\ \ (R_{\rm jup})$"]
label = [r"$100\times\,{\rm solar}$",
         r"$1.0\times\,{\rm solar}$",
         r"$0.1\times\,{\rm solar}$"]

# Figure 5:
h.histogram(Zp, Zr, parname=pname, percentile=[0.9545, 0.683], xran=xran,
   label=label, savefile="../plots/marginal-vs-metallicity.ps")


