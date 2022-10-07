import sys, glob
sys.path.append("/home/tb/secruz/DATA/PhD/P1_all_atom/step2/nadine-data/examples/TREMD/TREMD_kinetics/bin/kinetics")
import ala_kinetics
from ala_kinetics import *

import seaborn as sns
import sys, glob
import pandas as pd


name=sys.argv[1]

r=pd.read_pickle(name+"_100k_rates.pickle")
p=pd.read_pickle(name+"_100k_pop.pickle")
ev=pd.read_pickle(name+"_100kevents.pickle")
tp=pd.read_pickle(name+"_100k_tp.pickle")
tba=pd.read_pickle(name+"_100k_tba.pickle")
s = [[0.0, 1.0]]
rs = sym_counts_calc_rate(r, p, s, r.type.unique(), 100000)
#, time_unit_factor=1000.0)
r_s_ln = err_log_rate(rs, weight_name="sym_weight", diff_from_est=True )

unbinding = r[r.type == (1,0)][['temperature', 'rate', 'events', 'sum_weight']]
binding = r[r.type == (0,1)][['temperature', 'rate',  'events', 'sum_weight']]

np.savetxt('RESULTS/rates_unbinding.dat', unbinding.values)
np.savetxt('RESULTS/rates_binding.dat', binding.values)
np.savetxt('RESULTS/populations.dat', p.values)

error_unbinding = r_s_ln[r_s_ln.type == (1,0)][['temperature', 'rate', 'std_p', 'std_m', 'err_m', 'err_p' ]]
error_binding = r_s_ln[r_s_ln.type == (0,1)][['temperature', 'rate', 'std_p', 'std_m', 'err_m', 'err_p' ]]

np.savetxt('RESULTS/s-rates_unbinding.dat', error_unbinding.values)
np.savetxt('RESULTS/s-rates_binding.dat', error_binding.values)

tp_unbinding=tp[tp.type == (1,0)][['temperature', 'traj', 'start', 'stop', 'fraction']] 
tp_binding=tp[tp.type == (0,1)][['temperature', 'traj', 'start','stop','fraction']]

np.savetxt('RESULTS/tp_unbinding.dat', tp_unbinding.values)
np.savetxt('RESULTS/tp_binding.dat', tp_binding.values)
