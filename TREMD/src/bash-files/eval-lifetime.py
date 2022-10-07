import seaborn as sns
import sys, glob
sys.path.append("/home/tb/secruz/DATA/PhD/P1_all_atom/step2/nadine-data/examples/TREMD/TREMD_kinetics/bin/kinetics")
import ala_kinetics
from ala_kinetics import *
import sorted_lifetimes
from sorted_lifetimes import *

import pandas as pd

name=sys.argv[1]

r=pd.read_pickle(name+"_100k_rates.pickle")
p=pd.read_pickle(name+"_100k_pop.pickle")
ev=pd.read_pickle(name+"_100kevents.pickle")
tp=pd.read_pickle(name+"_100k_tp.pickle")
tba=pd.read_pickle(name+"_100k_tba.pickle")
s = [[0.0, 1.0]]


#calculate lifetimes for 0->1 (unbinding) transition & 1->0 (binding) transition
dw_01 = loop_dwell_trans_temp(tp, tba, [(0,1)])
dw_10 = loop_dwell_trans_temp(tp, tba, [(1,0)])

# Now adjust times according to weight for temperature 

dw_01_rw = dw_01[dw_01.temperature == 20 ].wait_T / dw_01[dw_01.temperature==20].weight
dw_10_rw = dw_10[dw_10.temperature == 20 ].wait_T / dw_10[dw_10.temperature==20].weight

np.savetxt('RESULTS/lifetime-binding.dat', dw_01_rw.values)
np.savetxt('RESULTS/lifetime-unbinding.dat', dw_10_rw.values)
