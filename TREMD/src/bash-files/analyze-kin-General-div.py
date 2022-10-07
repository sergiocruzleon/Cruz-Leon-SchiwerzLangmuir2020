import sys, glob
sys.path.append("/home/tb/secruz/DATA/PhD/P1_all_atom/step2/nadine-data/examples/TREMD/TREMD_kinetics/bin/kinetics/")
import ala_kinetics
from ala_kinetics import *

# Usage: python analyze-kin-General.py  ${ion} ${atom} ${path} ${temp_replica}
# Example python analyze-kin-General.py  Ba2 OP (O6G or N7) '/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/Ba2/rates-all/order-parameter/'  '/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/bash-files/index-temp/replica_temp_divalent.dat'
name=sys.argv[1]
atom=sys.argv[2]
path=sys.argv[3]
temp_repl=sys.argv[4]
max_num=sys.argv[5]
q_dict = {}
for i in range(0,20):
    fn = path+'/'+name+'/'+atom+'/'+atom+'_{}.xvg'.format(i)
    q_dict[i] = np.genfromtxt(fn, skip_header=0, comments='@' )[:int(max_num)]

#print first 5 lines from demux_0.xvg
print q_dict[0][:5]

#make np array from replica temp index
rep_ar = np.genfromtxt(temp_repl)[:int(max_num)]


#definition OP boundaries
s = [[0.0, 1.0]]

run_name=name+"_100k"

o_remd_st1 = get_rates_pop_multi_temp(q_dict, rep_ar, run_name, state_def=s,
                                      return_raw_events=True, remd=True,
                                      recalc=True, split_state_tp=True)

state_df_st1, all_tp_temp_df_st1, pt_st1, rates_st1, events_st1 = o_remd_st1


print rates_st1[rates_st1.temperature == 0]

print pt_st1[pt_st1.temperature == 0]

print all_tp_temp_df_st1[all_tp_temp_df_st1.temperature==0]

print rates_st1
