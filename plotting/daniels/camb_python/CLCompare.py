# sample script to compare two standard _lensedCl results
# ncl=1 does just TT, ncl=3,4 does more
# see camPlots.py for function declaration

from cambPlots import *
from pylab import *


root = r'../SPT+WMAP9_lensedCls.dat'
root2 = r'../Planck+WP_lensedCls.dat'


plot_compare([root, root2], ncl=1, x='sqrt')


show()

