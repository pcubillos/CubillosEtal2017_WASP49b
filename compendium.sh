# From the directory where this file is located, execute (this will make things easier):
topdir=`pwd`

# Clone (download) the PB code:
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 5d9e2f4

# Compile the PB code:
cd $topdir/pyratbay
make


# Download HITRAN/HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_w49b.txt
unzip '*.zip'
rm -f *.zip

# Make HITEMP H2O TLI file:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c tli_H2O.cfg


# Make atmospheric files:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c atm_000.10xsolar.cfg
$topdir/pyratbay/pbay.py -c atm_001.00xsolar.cfg
$topdir/pyratbay/pbay.py -c atm_100.00xsolar.cfg


# Make H2O opacity file:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c opacity_H2O.cfg


# Run MCMC for solar abundance model:
cd $topdir/run02/
$topdir/pyratbay/pbay.py -c mcmc_000.1xsolar.cfg
$topdir/pyratbay/pbay.py -c mcmc_001.0xsolar.cfg
$topdir/pyratbay/pbay.py -c mcmc_100.0xsolar.cfg

# High-resolution atmospheric runs:
cd $topdir/run03/
$topdir/pyratbay/pbay.py -c upper_atm_000.1x.cfg
$topdir/pyratbay/pbay.py -c upper_atm_001.0x.cfg
$topdir/pyratbay/pbay.py -c upper_atm_100.0x.cfg


# Figure 1:
cd $topdir
$topdir/fig_FORS2.py

# Figure 2:
cd $topdir
$topdir/fig_hydrodynamic.py

# Figure 3:
cd $topdir
$topdir/fig_radiative-eq.py

# Figure 4:
cd $topdir/run01/
$topdir/fig_clearspectra.py

# Figures 5 and 6:
cd $topdir/run02/
$topdir/fig_retrieval.py

# Figure 7:
cd $topdir
$topdir/fig_eqclouds.py

# Figure 8:
cd $topdir/run03/
$topdir/fig_hires.py
