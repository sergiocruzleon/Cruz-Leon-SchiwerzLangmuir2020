from ala_kinetics import *
from scipy.optimize import curve_fit

cl = sns.color_palette()
    
def dw_time_all_traj_w_ev(all_tp_temp, rep, first_time=0,
                          transition=None, return_prev_type=False):
    '''                      
    Setting transition to (0,0,0,1), enables me to look at the dwell time for 
    this transition. 
    '''
    col = ['temperature', 'type', 'traj', 'start', 'stop',
           'weight','wait_T', 'wait', 'prev_trans']
    df = pd.DataFrame(data=np.zeros((0, len(col))), columns=col)
    
    for trj_i in np.unique(all_tp_temp.traj.values):
        print trj_i
        df = dw_time_traj_w_ev(df, all_tp_temp[all_tp_temp.traj == trj_i], rep,
                               first_time=first_time, transition=transition,
                               return_prev_type=return_prev_type )
    return df


def dw_time_traj_w_ev(df, traj_all_tp_temp, rep, first_time=0, transition=None,
                      return_prev_type=False):
    '''
    The default for transition=None assumes that we have a two state system. 
    In a multi-state system we have to explicitly account for the fact that given state
    may be left in different ways. C.f. for instance Zwanzig 1983. 
    '''

    trj_i = traj_all_tp_temp.traj.unique()[0].astype(int) 
      
    #if transition:
    #   _sel_trans = traj_all_tp_temp[traj_all_tp_temp.type==transition]
    #   starts = _sel_trans.start.values
    #else:   
    starts = traj_all_tp_temp.start.values  
    stop = traj_all_tp_temp.stop.values

    start_u = np.unique(starts) # think about this
    stop_u = np.unique(stop)

    for s_i, start_tp in enumerate(start_u):
    
        # determine the type of the transition
        _tp = traj_all_tp_temp[traj_all_tp_temp.start ==start_tp]
        _tp_type = np.unique(_tp.type)   
        # Which temperatures contributed to the transition?

        
        if transition:
           #print _tp_type[0]
           if _tp_type[0] != transition:
              continue    
           
        _tp_temperatures =  np.unique(_tp.temperature.values)
          
    
        if s_i==0:
            end_preceding_tp = first_time  #q_r2_d[0][0,0]
        else:
             end_preceding_tp = stop_u[s_i -1]
        # I exclude frames that are already on the path 
        # add +1 to start, -1 to end_preceding_tp
    
        end_dw = start_tp -1
        start_dw = end_preceding_tp +1
        total_wait = end_dw - start_dw
        
        if return_prev_type and s_i > 0:
            _p = traj_all_tp_temp[traj_all_tp_temp.stop == start_dw-1]
            _p_type = np.unique(_p.type)[0][2:]
            #print _p_type
        else:
             _p_type = np.nan
             #print _p_type
    
        _temp_data = transition_path_data_pd(rep, start_dw, end_dw ).values[:, trj_i +1]
        frac = temp_frac_tp(_temp_data)
                
        for ti, w in frac.items():
        # exclude temperatures not sampled during the TP
            if ti in _tp_temperatures:
                t_weight = _tp[_tp.temperature == ti].fraction.values
        
                df = df.append({'temperature': ti, 'type': _tp_type[0], 'traj' : trj_i, 'start': start_dw,
                                'stop': end_dw, 'weight' : t_weight[0] , 'wait_T' : w* total_wait,
                                'wait': total_wait, 'prev_trans' : _p_type,}, ignore_index=True)
    return df
    

def get_bin_centres(b):
    return (b[:-1] + b[1:]) *0.5
    
    
def cumulative_hist(data, bins=100, weights=None):
    #print weights
    if weights is None:
        h, b = np.histogram(data, bins=bins,
                   density=False)
    else:
        h, b = np.histogram(data, bins=bins, weights=weights,
                           density=False)

    bc = get_bin_centres(b)
    cum = np.cumsum(h)
    cum_n = cum /  np.float64(cum.max())
    
    return bc, cum_n
   

def exp_dist_func(x, l):
    return l*np.exp(-l*x)


def cumulative_dist_func(x,l):
    return 1.0 - np.exp(-l*x)
    

def biexp_cdf(t, w, tau1, tau2):
    return 1 - w*np.exp(-t/tau1) - (1-w)*np.exp(-t/tau2)
    

def fit_waiting_cdf(data, bins=100, weights=None,
                   return_noncumulative=False, return_bin_edges=False,
                   fit_cumulative_cdf=True):
    """
    fit_func: cumlative_dist_func or exp_dist_func
    """
    
    bc, cum_n  = cumulative_hist(data, bins=bins, weights=weights)
    
    if fit_cumulative_cdf:
       popt = curve_fit(cumulative_dist_func, bc, cum_n, p0=0.002)
    else:
         popt = curve_fit(exp_dist_func, bc, cum_n, p0=0.002)   
    
    return_l = [bc, cum_n, popt[0]]
    
    if return_noncumulative:
        return_l.append(h)
    if return_bin_edges:
        return_l.append(b)
    
    return return_l

    
def fit_plot_cdf(ax, wait_T, dist_label='CDF', fit_label='fit CDF', bins=1000, weights=None,
                 return_h=False, plot_fit=True, color=None):
                 
    b, h, l = fit_waiting_cdf(wait_T, bins=bins, weights=weights)
    if color:
       ax.plot(b, h, lw=2, label=dist_label, c=color)
    else:   
        ax.plot(b, h, lw=2, label=dist_label)
    if plot_fit:
       ax.plot(b, cumulative_dist_func(b, l))

    if return_h:
       return ax, b, l, h
    else:   
         return ax, b, l


def format_cdf(fig, ax, time_factor=1.0, time_units_label="x 1000 MD steps"):
    ax.set_xlabel('Time '+time_units_label, fontsize=14)
    ax.set_xlim(10**1*time_factor, 10**5.5*time_factor)
    ax.set_ylabel("CDF", fontsize=14)
    fig.tight_layout()
    return fig, ax


def plot_cdf_multi_temp(dw, ax, temp_l, bins=100, temp_label_l=None, dist_label_l=None,
                        rate_type=None, rate_factor=1):
                        
    
    for ti, temp in enumerate(temp_l):
        _dw = dw[dw.temperature == temp]
        if temp_label_l:
           t_label=' T={:.0f} '.format(temp_label_l[ti])
        else:
             t_label=' T={:.0f}'.format(temp)
        if dist_label_l:
           dist_label = dist_label_l[ti] + t_label
        else:
             dist_label = t_label   
        ax[ti], b, l = fit_plot_cdf(ax[ti], _dw.wait_T / _dw.weight, bins=bins, dist_label=dist_label)
        ax[ti].legend(loc=8)
        
        if rate_type is not None:
           _r = rate_type[rate_type.temperature == temp].rate.values[0] / rate_factor
           ax[ti].plot(b, cumlative_dist_func(b, _r), '--') # , label='rate coeff'
    return ax


def plot_cdf_tp(dw, ax, temp_l, bins=100, temp_label_l=None, dist_label_l=None,
                rate_type=None, rate_factor=1, plot_fit=True):
    
    for ti, temp in enumerate(temp_l):
        _dw = dw[dw.temperature == temp]
        _tp = _dw.stop - _dw.start
        #print _tp.head()
        if temp_label_l:
           t_label=' T={:.0f} '.format(temp_label_l[ti])
        else:
             t_label=' T={:.0f}'.format(temp)
        if dist_label_l:
           dist_label = dist_label_l[ti] + t_label
        else:
             dist_label = t_label   
        ax[ti], b, l = fit_plot_cdf(ax[ti], _tp, weights=_dw.fraction.values,  bins=bins, dist_label=dist_label,
                                    plot_fit=plot_fit)
        ax[ti].legend(loc=8)
        
        if rate_type is not None:
           _r = rate_type[rate_type.temperature == temp].rate.values[0] / rate_factor
           ax[ti].plot(b, cumlative_dist_func(b, _r), '--')
    return ax


def _tp_length():
    pass     
# new attempt


def plot_compare_lifetime_distributions(ax, bins, dist1, dist2, dist1_label, dist2_label, dist1_cl, dist2_cl,
                                         hw = 2 * 10, hl = 0.04 * 2, arrow=None, arrow_sign=1, verbose=False,
                                         arrow_start=0.22, exclude_empty_waits=True):
    """

    :param ax:
    :param bins:
    :param dist1:
    :param dist2:
    :param dist1_label:
    :param dist2_label:
    :param dist1_cl:
    :param dist2_cl:
    :param _hw:
    :param _hl:
    :param arrow: None, mean, median, fit
    :param arrow_sign:
    :param verbose:
    :return:
    """
    ax.semilogx()

    ax, b_md, l_md, hist_md_dw = fit_plot_cdf(ax, dist1, plot_fit=False, color=dist1_cl,
                                             return_h=True, bins=bins, dist_label=dist1_label)

    ax, b_remd, l_remd, hist_remd_dw = fit_plot_cdf(ax, dist2, color=dist2_cl,
                                                   return_h=True, plot_fit=False, bins=bins,
                                                   dist_label=dist2_label)
    if exclude_empty_waits:
       dist1_ = dist1[dist1 > 0]
       dist2_ = dist2[dist2 > 0]
    else:
        dist1_ = dist1
        dist2_ = dist2
    
    if arrow == "mean":
        arrow_md = np.mean(dist1_)
        arrow_remd = np.mean(dist2_)
    elif arrow == "median":
        arrow_md = np.median(dist1_)
        arrow_remd = np.median(dist2_)
    elif arrow == "fit":
        arrow_md = l_md
        arrow_remd = l_remd

    if arrow:
        print (1/arrow_md, 1/arrow_remd)

    if  arrow == "fit":
        ax.arrow(1.0 / arrow_md, 1 - np.exp(-1), 0, arrow_sign*0.22, color=dist1_cl, head_width=hw, head_length=hl,
                 lw=2)
        ax.arrow(1.0 / arrow_remd, 1 - np.exp(-1), 0, arrow_sign*0.22, color=dist2_cl, head_width=hw, head_length=hl,
                  lw=2)
    
    elif arrow == "mean":
         ax.arrow(arrow_md, 1 - np.exp(-1), 0, arrow_sign*0.22, color=dist1_cl, head_width=hw, head_length=hl,
                 lw=2)
         ax.arrow(arrow_remd, 1 - np.exp(-1), 0, arrow_sign*0.22, color=dist2_cl, head_width=hw, head_length=hl,
                  lw=2)
    elif arrow == "median":
         ax.arrow(arrow_md, 0.5, 0, arrow_sign*0.22, color=dist1_cl, head_width=hw, head_length=hl, lw=2)
         ax.arrow(arrow_remd, 0.5, 0, arrow_sign*0.22, color=dist2_cl, head_width=hw, head_length=hl, lw=2)

    return ax
