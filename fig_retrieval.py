#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.plots as mp
import MCcubed.utils as mu

sys.path.append("../pyratbay")
import pyratbay.constants  as pc
import pyratbay.atmosphere as atm

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
def readmolfile(molfile):
  """
  Read Pyrat Bay's molecular info file.
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


def histogram(Z1, Z2=None, Z3=None, label=None, parname=None, thinning=1,
               savefile=None, percentile=[], xran=None):
  """
  Plot marginal posteriors for a set of parameters (columns) and runs (rows).
  """
  nhists = len(Z1)
  if Z2 is None:
    npost = 1
  elif Z3 is None:
    npost = 2
  else:
    npost = 3

  nh = npost*nhists
  pdf  = [None]*nh
  xpdf = [None]*nh
  if not isinstance(pdf, list):  # Put single arrays into list
    pdf  = [pdf]
    xpdf = [xpdf]
  # Histogram keywords depending whether one wants the HPD or not:
  hkw = {}
  nper = 0
  if percentile is not None:
    hkw = {'histtype':'step', 'lw':2, 'color':'b'}
    nper = len(percentile)
    boundaries = np.zeros((nhists, nper, 2), np.double)

  fs = 14  # Fontsize

  nrows    = nhists
  ncolumns = npost

  histheight = np.amin((2 + 2*(nrows), 8))
  plt.figure(-2, figsize=(8, histheight))
  plt.clf()
  plt.subplots_adjust(left=0.15, right=0.95, bottom=0.18, top=0.95,
                      hspace=0.07, wspace=0.1)
  post = [Z1, Z2, Z3]
  for j in np.arange(npost):
    maxylim = 0  # Max Y limit
    for i in np.arange(nhists):
      ax = plt.subplot(nhists, ncolumns, npost*i+j+1)
      a  = plt.xticks(size=fs-1.5, rotation=90)
      if j == 0:
        a = plt.yticks(size=fs-1.5)
      #else:
      a = plt.yticks(visible=False)
      if i+1 == nhists:
        plt.xlabel(parname[j], size=fs)
      else:
        a = plt.xticks(visible=False)
      vals, bins, h = plt.hist(post[j][i][0::thinning], bins=25,
                               normed=True, range=xran[j], **hkw)
      vals = np.r_[0, vals, 0]
      bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
      # Plot HPD region:
      k = 1
      fc = ['0.65', '0.9']
      for p in np.arange(nper):
      #if percentile is not None:
        PDF, Xpdf, HPDmin = mu.credregion(post[j][i][:], percentile[p],
                                          pdf[npost*i+j], xpdf[npost*i+j])
        # interpolate xpdf into the histogram:
        f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
        # Plot the HPD region as shaded areas:
        xr = np.argwhere((Xpdf>xran[j][0]) & (Xpdf<xran[j][1]))
        Xpdf = Xpdf[np.amin(xr):np.amax(xr)]
        PDF  = PDF [np.amin(xr):np.amax(xr)]
        # Boundaries of the HPD region:
        limXpdf = np.amin(Xpdf[PDF>HPDmin]), np.amax(Xpdf[PDF>HPDmin])
        if j == 0: # Only for transit-pressure posterior
          boundaries[i,p] = limXpdf
        ax.fill_between(Xpdf, 0, f(Xpdf),
          where=((Xpdf>limXpdf[0]) & (Xpdf<limXpdf[1])),
          facecolor=fc[nper-k], edgecolor='none', interpolate=False)
        k += 1

      ax.set_xlim(xran[j][0], xran[j][1])
      maxylim = np.amax((maxylim, ax.get_ylim()[1]))

    # Set uniform height:
    maxylim *= 1.1
    for i in np.arange(nhists):
      ax = plt.subplot(nrows, ncolumns, npost*i+j+1)
      ax.set_ylim(0, maxylim)
      if j == 0:
        plt.text(-1.3, 0.85*maxylim, label[i], ha="right", fontsize=15)
      elif j==1:
        plt.plot([1.198, 1.198], [0,maxylim], "r--", lw=2.5)

  if savefile is not None:
    plt.savefig(savefile)
  # Save HPD boundaries to file:
  with open("pt_boundaries.dat", "w") as f:
    for p in np.arange(nper)[::-1]:
      f.write("{:.1f}%          ".format(100.0*percentile[p]))
    f.write("\n")
    for i in np.arange(nhists):
      for p in np.arange(nper)[::-1]:
        f.write("{:5.2f} {:5.2f}    ".format(boundaries[i,p,0],
                                             boundaries[i,p,1]))
      f.write("\n")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ID, symbol, mass, diam = readmolfile("../pyratbay/inputs/molecules.dat")

p0       = 0.1   * pc.bar   # Reference pressure level
rtransit = 1.198 * pc.rjup  # Transit radius (white light)
mplanet  = 0.396 * pc.mjup  # Planetary mass
rstar    = 1.038 * pc.rsun  # Stellar radius

atmfile = ["../run01/WASP-49b_100.0xsolar.atm",
           "../run01/WASP-49b_001.0xsolar.atm",
           "../run01/WASP-49b_000.1xsolar.atm"]

postfiles = ["./MCMC_WASP-49b_100.0xsolar_Tprior.npz",
             "./MCMC_WASP-49b_001.0xsolar_Tprior.npz",
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
  mm      = np.sum(q*mspec, axis=1) # Mean molecular mass profile
  
  # Unique posterior values:
  utemp, uind, uinv = np.unique(temp, return_index=True, return_inverse=True)
  urad   = rad[uind]
  upress = np.zeros(len(utemp), np.double)
  
  nsmall = 0
  for i in np.arange(len(utemp)):
    # Compute radius atmospheric profile:
    radius = atm.hydro_m(press, utemp[i], mm, mplanet, p0, urad[i])
    if radius[0] < rtransit:
      nsmall += 1
      continue
    # Interpolate pressure-radius:
    f = si.interp1d(radius, logp)
    # Get pressure of the photospheric radius (transit radius):
    upress[i] = 10**f(rtransit) / pc.bar
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
  # Pair-wise plots:
  mp.pairwise(Z, parname=parname, fignum=j+10,
              savefile="../plots/pair_{}.ps".format(name[j]))
  mp.histogram(Z, parname=parname, percentile=0.683, fignum=20+j)

xran = [[-8, -1.0], [1.06, 1.205]]
pname = [r"$\log_{10}(p_{\rm T})\ \ ({\rm bar})$",
         r"$R_{0.1\,{\rm bar}}\ \ (R_{\rm jup})$"]
label = [r"$100\times\,{\rm solar}$",
         r"$1.0\times\,{\rm solar}$",
         r"$0.1\times\,{\rm solar}$"]

# Compbined marginal posterior histograms plot:
boundaries = histogram(Zp, Zr, parname=pname, percentile=[0.9545, 0.683],
            xran=xran, label=label, savefile="../plots/marginal_posteriors.ps")


