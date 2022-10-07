import numpy as np
import collections
from collections import Counter
import glob
import os
from copy import deepcopy
from collections import defaultdict

def get_rates_pop(q_ar, f_centre=0.6, uf_centre=0.1):
    # define starting point
    
    f_events = []
    uf_events = []
    
    out_ar = np.zeros(q_ar.shape)
    out_ar[:,0] = q_ar[:,0]
   # print out_ar, out_ar.shape
    
    if np.abs(q_ar[0,1] - uf_centre) < np.abs(q_ar[0,1] - f_centre):
        out_ar[0,1] =  0.1 #0.0
    else:
         out_ar[0,1] = 0.6 #1.0
    
    for i, row in enumerate(q_ar[1:,:]):
        if row[1] >= f_centre:
            v = 0.6 # 1.0
        elif row[1] <= uf_centre:
            v = 0.1 #0.0
        else: # v stays the same
            v = out_ar[i, 1]
            
        out_ar[i+1,1] = v
        
        assert( q_ar[i+1,0] == out_ar[i+1,0]),  "{} {}".format(q_ar[i,0],out_ar[i+1,0])
        
        if v > out_ar[i,1]:
            f_events.append(row[0])
        elif v < out_ar[i,1]:
            uf_events.append(row[0])
    
    return out_ar, f_events, uf_events

def find_transition_paths(q_ar, f_centre=0.6, uf_centre=0.1):
    """

    :param q_ar:
    :param f_centre:
    :param uf_centre:
    :return: out_ar, f_events, uf_events

    Finds continuous paths between f_centre (folded or native state) and uf_centre (unfolded state), which
    are the transitions paths. The function returns an array with each time point assigned to either native or
    the unfolded state. The function also returns a list
    """

    f_events = []
    uf_events = []
    tps_start_time = 0.0
    on_tps = False

    out_ar = np.zeros(q_ar.shape)
    out_ar[:,0] = q_ar[:,0]
   # print out_ar, out_ar.shape

    if np.abs(q_ar[0,1] - uf_centre) < np.abs(q_ar[0,1] - f_centre):
        out_ar[0,1] =  0.1 #0.0
    else:
         out_ar[0,1] = 0.6 #1.0

    for i, row in enumerate(q_ar[1:,:]):
        if row[1] >= f_centre:
            v =0.6
            on_tps = False
        elif row[1] <= uf_centre:
            v = 0.1
            on_tps = False
        else:
            v = out_ar[i,1]
            # if not in one of the two endstates, the molecule may be on a transition
            if not on_tps:
                tps_start_time = row[0] # save the time in case a true transition happens
                on_tps = True

        out_ar[i+1,1] = v

        if v > out_ar[i,1]:
            assert (on_tps == False)
            f_events.append((tps_start_time, row[0]))
        elif v < out_ar[i,1]:
            assert (on_tps == False)
            uf_events.append((tps_start_time, row[0]))
    return out_ar, f_events, uf_events

def temp_transition_path(rep_ar, event_l, trj_index, verbose=False):
    """

    :param rep_ar:
    :param event_l:
    :param trj_index:
    :param verbose:
    :return:
    """
    transition_temperature_d = {}

    #for ev in event_l
    # events list is ordered

    events_counter=0
    on_tp=False
    temperatures_tp = []
    times_tp = []

    #for row in rep_ar[1:,:]:
    for row in rep_ar:
        if on_tp:
            if row[0] == event_l[events_counter][1]:
                on_tp =False
                temperatures_tp.append(row[trj_index + 1])
                times_tp.append(row[0])
                transition_temperature_d[(trj_index, event_l[events_counter])] = np.column_stack((times_tp,
                                                                                                  temperatures_tp))

                if verbose:
                    print event_l[events_counter]
                    print temperatures_tp

                events_counter += 1

                if events_counter == len(event_l):
                    break
            else:
                 temperatures_tp.append(row[trj_index + 1])
                 times_tp.append(row[0])
        #print row[0], events_counter

        elif row[0] == event_l[events_counter][0]:
            temperatures_tp = [row[trj_index + 1]]
            times_tp = [row[0]]
            on_tp = True
    if verbose:
       print events_counter, len(event_l) # assert this?
    return transition_temperature_d

def tp_length(events_l):
    tp_time_l = []
    for ev in events_l:
        #print ev[1] - ev[0]
        tp_time_l.append(ev[1] - ev[0])
    return  tp_time_l

# the output_array should have the same shape as the replica index array, which should enable me to assign the time
# spent in folded and unfolded states at each temperature

def event_prob(events_d, number_traj): # this does not work
    """
    :param events_d:
    :param number_traj:
    :return:

    prob_d is ordered by trajectories not by temperatures?

    """
    prob_d = {}
    for n in range(number_traj):
        #print k
        prob_d[int(n)] = []

    for k, v in events_d.items():
        tp_length = len(v[:, 1])
        #print v[1]
        for temp, c in Counter(v[:,1]).items():
        # http://stackoverflow.com/questions/2600191/how-can-i-count-the-occurrences-of-a-list-item-in-python
           # print temp, c
            prob_d[int(temp)].append(c / np.float64(tp_length))
    return prob_d

def sum_temp_prob_demux_trj(prob_d):
    demux_tp_temperature = {}
    for k, v in prob_d.items():
        demux_tp_temperature[k] = np.sum(v)
    return demux_tp_temperature

def sum_temp_prob_demux_trj_ar(prob_d, trj_ar, trj_i):
    '''
    Transition probabilities for each temperature for each trajectory

    each row is one trajectory, each column one temperature

    '''
    for k, v in prob_d.items():
        trj_ar[trj_i, k] = np.sum(v)
    return trj_ar

def event_prob_demux(events_d, number_temperatures):
    prob_d = event_prob(events_d, number_temperatures)
    return sum_temp_prob_demux_trj(prob_d)

def analyse_demux_trj(q_ar, repl_ar, trj_i, f_ar, uf_ar, tp_tuple=None):
    """

    Detect all transition events for each trajectory.

    Loop through the transition events and assign fractional transition probabilities for each temperature.

    :param q_ar:
    :param repl_ar:
    :param trj_i:
    :param f_ar:
    :param uf_ar:
    :return:
    """
    out_ar, f_events, uf_events = find_transition_paths(q_ar)

    if tp_tuple:
        tp_reaction_coord_d = tp_tuple[0]
        tp_temperature_d = tp_tuple[1]
        folding_unfolding_event_d = tp_tuple[2]

    for i, event in enumerate([[f_events, f_ar], [uf_events, uf_ar]]):

        _tp_rc = extract_transition_path_data(event[0], q_ar, trj_i)

        event_d = temp_transition_path(repl_ar, event[0], trj_i)
        event_temp_prob = event_prob(event_d, len(repl_ar[0,1:]))
        event[1] = sum_temp_prob_demux_trj_ar(event_temp_prob, event[1], trj_i) # changing a array via its list entry

        if tp_tuple:
            #print trj_i
            tp_reaction_coord_d.update(_tp_rc)
            tp_temperature_d.update(event_d)
            # order event probability temperature dictionaries by temperature indices
            folding_unfolding_event_d[trj_i, i] = event[0]

    if tp_tuple:
        return out_ar[:,1], f_ar, uf_ar, (tp_reaction_coord_d, tp_temperature_d, folding_unfolding_event_d)
    else:
        return out_ar[:,1], f_ar, uf_ar


    # prob_folding_temp = sum_temp_prob_demux_trj(f_d)
   # prob_unfolding_temp =sum_temp_prob_demux_trj(uf_d)
   #return sum_temp_prob_demux_trj_ar(prob_d, trj_ar, trj_i), sum_temp_prob_demux_trj_ar(prob_d, trj_ar, trj_i)

# function to run the analysis of REMD; also add function to analyse straightforward MD

def kinetics_from_remd(q_dict, rep_ar, dt=5.0, output_folding_events=False, output_transition_paths=False,
                       start_time=0):
    """

    :param q_dict:
    :param rep_ar:
    :param dt:
    :param output_out_ar:
    :param output_folding_events:
    :return:

    The function loops through all REMD trajectories, in which temperature is a dynamic variable. For each trajectory
    the transition events are found. The event probabilities are distributed across the temperatures at which the event
    occurred.
    """
    if output_transition_paths:
        tp_reaction_coord_d = {}
        tp_temperature_d = {}
        folding_unfolding_event_d = {}
        tp_tuple = (tp_reaction_coord_d, tp_temperature_d, folding_unfolding_event_d)
    else:
        tp_tuple = None

    num_temp = rep_ar[0, 1:].size
    rep_ar = rep_ar[rep_ar[:,0] >= q_dict[0][start_time,0]]
    num_frames = len(rep_ar)
    out_ar = np.zeros((num_frames, num_temp))

    folding_events_p_ar = np.zeros((num_temp, num_temp), dtype=np.float64)
    unfolding_events_p_ar = np.zeros((num_temp, num_temp), dtype=np.float64)

    # put into function that can return different things?
    for trj_i, q_fn in q_dict.items():
        q_ar = q_fn
        q_ar = q_ar[start_time:, :]
        _demux = analyse_demux_trj(q_ar, rep_ar[start_time:, :], trj_i,
                                   folding_events_p_ar, unfolding_events_p_ar, tp_tuple=tp_tuple)

        if output_transition_paths:
            out_ar[:, trj_i], folding_events_p_ar, unfolding_events_p_ar, tp_tuple = _demux
        else:
            out_ar[:, trj_i], folding_events_p_ar, unfolding_events_p_ar = _demux

   # out_ar_temp = reorder_pop_temp(rep_ar, out_ar)
   # pf_t, pu_t = time_state(out_ar_temp, 0.6), time_state(out_ar_temp, 0.1)

    out_ar_temp = reorder_temperature_indices(out_ar, rep_ar, arrange_up_to=num_temp)
    out_ar_bin = exchange2binary_val(out_ar_temp)
    pf_t = np.mean(out_ar_bin, axis=0)
    pu_t = 1.0 - pf_t
    rate_folding_t = (np.sum(folding_events_p_ar, axis=0)) / (pu_t * rep_ar[-1,0])
    rate_unfolding_t = (np.sum(unfolding_events_p_ar, axis=0)) / (pf_t * rep_ar[-1,0])

   # return pf_t/np.float64(num_frames), pu_t/np.float64(num_frames), rate_folding_t, rate_unfolding_t

    if output_folding_events:
        if output_transition_paths:
            return out_ar_bin, pf_t, pu_t, rate_folding_t, rate_unfolding_t,\
                   [folding_events_p_ar, unfolding_events_p_ar], tp_tuple
        else:
            return out_ar_bin, pf_t, pu_t, rate_folding_t, rate_unfolding_t,\
                   [folding_events_p_ar, unfolding_events_p_ar]
    elif output_transition_paths:
        return out_ar_bin, pf_t, pu_t, rate_folding_t, rate_unfolding_t, tp_tuple
    else:
        return out_ar_bin, pf_t, pu_t, rate_folding_t, rate_unfolding_t

def calc_rates_events(events, p_react, length_trj):
    """
    :param events:
    :param p_react:
    :param length_trj:
    :return:
    """
    return events / (p_react * length_trj)

def reorder_pop_temp(rep_ar, out_ar):
    """
    :param rep_ar:
    :param out_ar:
    :return:

    NB needs a unit test
    DOES NOT SEEM TO WORK/ could a problem of usage
    """
    rep_int = np.zeros(rep_ar[:, 1:].shape)
    rep_int[:] = rep_ar[:, 1:] #.astype(np.int)
    rep_int = rep_int.astype(np.int)
    out_ar_temp = np.zeros(out_ar.shape)

    for ri, row in enumerate(out_ar):
        out_ar_temp[ri, :] = row[rep_int[ri]]
    return out_ar_temp

def reorder_temperature_indices(out_ar, rep_ar, arrange_up_to=24):
    """

    :param out_ar:
    :param rep_ar:
    :return:

    The function loops through the temperatures, All frames (i.e. reaction coordinate values) for a temperature are
    written into a column. The selection of the temperatures based on the indices in the rep_ar has been tested for
    a minimal example, which needs to be formalised into a unit test.

    This seems to work. Need a unit test
    """
    test_reorder = np.zeros(out_ar.shape)
    ti_l = np.arange(0, arrange_up_to)
    for ti in ti_l:
        test_reorder[:, ti] = out_ar[:, :][rep_ar[:, 1:] == ti]
    return test_reorder


def time_state(data_ar, cut_off):
    return np.sum(data_ar == cut_off, axis=0)

def extract_transition_path_data(event_l, data_ar, trj_index, verbose=False):
    """
    Extract the value of the reaction coordinate during a transition P(Q| TP)

    :param event_l:
    :param data_ar:
    :param trj_index:
    :param verbose:
    :return:
    """

    transition_path_data_d = {}
    events_counter = 0
    on_tp = False
    tp_data_l = []

    for row in data_ar:
        if on_tp:
            if row[0] == event_l[events_counter][1]:
                on_tp = False
                tp_data_l.append(row)
                transition_path_data_d[(trj_index, event_l[events_counter])] = tp_data_l # events_counter

                if verbose:
                    print event_l[events_counter]
                    print tp_data_l

                events_counter += 1

                if events_counter == len(event_l):
                    break
            else:
                tp_data_l.append(row)

        elif row[0] == event_l[events_counter][0]:
            tp_data_l = [row]
            on_tp = True
    if verbose:
       print events_counter, len(event_l) # assert this?
    return transition_path_data_d


def make_q_dict(path, split_char, first_index=0, last_index=None, split_pos=-1):
    q_l = glob.glob(path)
    q_dict = {}
    for q_fn in q_l:
        trj = os.path.basename(q_fn).split(split_char)[split_pos] # "_c"
        q_dict[int(float(trj))] = np.genfromtxt(q_fn)[first_index:last_index,:]
    return q_dict


def exchange2binary_val(test_reorder,f_centre=0.6, uf_centre=0.1, new_folded_val=1.0, new_unfolded_val=0.0):
    """

    :param test_reorder:
    :param f_centre:
    :param uf_centre:
    :return:

    Exchange the value in the arrays with the state assignments for zeros and ones.
    """
    out_bin = deepcopy(test_reorder)
    out_bin[ out_bin == f_centre] = new_folded_val
    out_bin[out_bin == uf_centre] = new_unfolded_val
    return out_bin

def reorder_tp_by_temperature(tp_temp, tp_rc):
    tp_temp_q_dict = defaultdict(list)

    for k, v in tp_temp.items():

        for i, t in enumerate(v[:, 1]):
            _rc = tp_rc[k][i]
            tp_temp_q_dict[t].append(list(_rc))
    return tp_temp_q_dict


def get_rates_md(path, split_char, first_index=0, last_index=None, verbose=False):

    q_dict = make_q_dict(path, split_char, first_index=first_index, last_index=last_index )

    t_out_d = {}
    t_f_ev = {}
    t_u_ev = {}

    for qi, q_ar in q_dict.items():
       # print qi, q_ar
        out_ar, f_events, uf_events = find_transition_paths(q_ar)
        t_u_ev[qi] = uf_events
        t_f_ev[qi] = f_events
        t_out_d[qi] = out_ar

    #print q_dict
    sorted_keys = sorted(q_dict.keys())
    print sorted_keys

    temp_ind_d = {}
    for ki, key in enumerate(sorted_keys):
        temp_ind_d[key]= ki

    pf_temp_d = {}
    puf_temp_d = {}
    f_rate_d = {}
    uf_rate_d = {}
    #temp_pop_ar = np.zeros((len(q_dict.keys()), 2))
    for temp, out_ar in t_out_d.items():
        out_bin = exchange2binary_val(out_ar)
        p_fold_temp = np.sum(out_bin[:, 1], axis=0) / len(out_bin)
        #temp_pop_ar[temp_ind_d[temp], :] = temp, p_fold_temp
        pf_temp_d[temp] = p_fold_temp
        p_ufold_temp = 1.0 - p_fold_temp
        puf_temp_d[temp] = p_ufold_temp
        f_rate_d[temp] = calc_rates_events(len(t_f_ev[temp]), p_ufold_temp, out_ar[-1,0])
        uf_rate_d[temp] = calc_rates_events(len(t_u_ev[temp]),p_fold_temp, out_ar[-1,0])

    #folding_rate_t = calc_rates_from_events_md(t_f_ev, puf_temp_d)
    #unfolding_rate_t = calc_rates_from_events_md(t_u_ev, pf_temp_d)
    #return folding_rate_t, unfolding_rate_t
    return collections.OrderedDict(sorted(f_rate_d.items())),\
           collections.OrderedDict(sorted(uf_rate_d.items()))

# http://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key

def calc_rates_from_events_md(event_temp_d, pfold_temp):
    """

    :param event_temp_d:
    :param pfold_temp:
    :return:
    """
    #temp_rates_ar = np.zeros((len(event_temp_d.keys()),2))
    temp_rates_d = {}
    for t, events_l in event_temp_d.items():
        print
        traj_length =[-1,0]
        temp_rates_d[t] = t, calc_rates_events(len(events_l), pfold_temp[t][:, 1], traj_length)
    return collections.OrderedDict(sorted(temp_rates_d.items()))

