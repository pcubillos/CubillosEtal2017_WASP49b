# Clone code
# From the directory where this file is located, execute:
topdir=`pwd`
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 0b47988
make


# TEA patch: with a text editor open:
$topdir/pyratbay/modules/TEA/tea/balance.py
# And change line 147 to:
            free_id.append(n + m)


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

# Figure 1:
cd $topdir
$topdir/fig_FORS2.py

# Figure 4:
cd $topdir/run01/
$topdir/fig_clearspectra.py

# Figures 5 and 6:
cd $topdir/run02/
$topdir/fig_retrieval.py


