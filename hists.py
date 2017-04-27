import sys, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as si

here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(here + '/pyratbay/modules/MCcubed')
import MCcubed.utils as mu


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

  with open("pt_boundaries.dat", "w") as f:
    for p in np.arange(nper)[::-1]:
      f.write("{:.1f}%          ".format(100.0*percentile[p]))
    f.write("\n")
    for i in np.arange(nhists):
      for p in np.arange(nper)[::-1]:
        f.write("{:5.2f} {:5.2f}    ".format(boundaries[i,p,0],
                                             boundaries[i,p,1]))
      f.write("\n")
