__author__ = 'lustelzl'
from kinetics.ala_kinetics import *
from numpy.testing import assert_array_equal
from kinetics.sorted_lifetimes import *
from hummer import transition
from copy import deepcopy

class TestKineticsAla:

    def __init__(self):
        self.ar = np.array([[0.0,  0.3, 0, 0.09, 0.1, 0.09, 0.11, 1.0, 1.0],
                            [0.4, 0.5, 1, 1.00, 0.0, 0.0, 0.1, 1.0, 1.0]])

        self.t_pts = np.arange(0,len(self.ar[0, :]))
        self.r = np.column_stack((self.t_pts.T, self.ar.T))
        self.s = [[0.1, 0.6], [0.1, 0.6]]

        self.trj2 = np.zeros(self.r.shape)
        self.trj2[:, 1:] = 0.01
        self.trj2[:,0] = self.t_pts
        self.trj2[1:5, 1] = [0.65] * self.r[1:5, 1].size
        self.trj2[5:, 2] = [0.65] * self.r[5:, 1].size
        
        self.remd, self.rep_ar, self.sdf = self.setup_remd_trj()
        self.mod_remd = self.setup_mod_remd_trj()
        
        o = get_rates_pop_multi_temp(self.remd, self.rep_ar, "test_remd",
                                     state_def=self.sdf, return_raw_events=True)
        self.state_df, self.all_tp_temp_df, self.pt, self.rates, self.events_d  =o


    def test_tba_2coord(self):
        self.setup_tba(self.r, self.s)
        self.setup_tba(self.trj2, self.s)


    def setup_tba(self, ar, borders):
        t_bin, event_d = multi_tp(ar, borders, verbose=True)

        tba1 = transition.tba(ar[:, 1], borders[0][0], borders[0][1])
        tba2 = transition.tba(ar[:, 2], borders[1][0], borders[1][1])
        assert_array_equal(t_bin[:, 1], tba1)
        assert_array_equal(t_bin[:, 2], tba2)

# need to test event time identification

# test calcuation of rates from events at multiple temperatures

# test calculation of block averages

    def setup_remd_trj(self):
        """
        hard coded REMD example trajectory for testing 
        """
        r_t = np.arange(10)
        trj1 = np.array([[0,0,],
                        [0,0.11],
                        [0,0.2],
                        [0,0.7],
                        [0,0.8]])
        trj1 = np.vstack((trj1, trj1[::-1]))
        trj1 = np.column_stack(( r_t, trj1,))
        
        trj2 = np.array([[0.8, 0.7],
                        [0.05, 0.02],
                        [0.05, 0.02],
                        [0.0, 0.0],
                        [0.1, 0.1],
                        [0.65, 0.8],
                        [0.8, 0.65],
                        [0.05, 0.02],
                        [0, 0],
                        [0, 0]])
        trj2 = np.column_stack(( r_t, trj2))
        
        t = np.array([[1, 0], [1,0], [0,1], [1,0], [0,1]])
        t_t = np.column_stack((np.arange(10), np.vstack((t, t[::-1]))))
        s= [[0.1, 0.6], [0.1,0.6]]
        trj_dict = {}
        trj_dict[0] = trj1
        trj_dict[1] = trj2
        
        return trj_dict, t_t, s
    
    
    def setup_mod_remd_trj(self):
        mod_remd = deepcopy(self.remd)
        mod_remd[0][3, -1] = 0.59
        mod_remd[0][8, -1] = 0.0  
        return mod_remd 
       
    def _test_plot(self, state_df, all_tp_temp_df, save_fig=True):
        """Plotting the test case""" 
        
        import matplotlib as mpl
        cl = sns.color_palette("coolwarm", 2)
        cl_cmap = mpl.colors.ListedColormap(cl)
        
        #state_df, all_tp_temp_df = run_multi_tp(self.remd, self.rep_ar,self.sdf)
        
        s_trj1 = state_df.loc[0][['a', 'b']].apply(state_bin2num, axis=1)
        s_trj2 = state_df.loc[1][['a', 'b']].apply(state_bin2num, axis=1)
        
        trj1 = self.remd[0]
        trj2 = self.remd[1]
        
        fig, ax = plt.subplots(1,2, figsize=(10,3.5))

        ax[0].plot(trj1[:,0], trj1[:,1], ".")
        ax[0].plot(trj1[:,0], trj1[:,2], ".")

        ax[0].set_ylim(-0.1, 1.1)
        #ax[0].set_xlim(-0.2, 3.2)

        ax[0].plot([0,10], [0.1]*2, "--")
        ax[0].plot([0,10], [0.6]*2, "--")

        #ax[0].set_xticks(self.rep_ar[0,1:])
        #ax[1].set_xticks(self.rep_ar[0,1:])
        ax[1].set_ylim(-0.2,3.2)

        ax[1].scatter(state_df.loc[0].index.values.astype(np.float64), s_trj1.values,
                      c= state_df.loc[0]['temperature'].values, cmap = cl_cmap, s=50 )

        ax[1].plot(s_trj1, "-", lw=0.5, c='grey')
        ax[1].set_yticks([0,1,2,3])
        
        if save_fig:
           fig.savefig("trj1_remd.png")
        
        # plotting trajectory 2
        fig, ax = plt.subplots(1,2, figsize=(10,3))

        ax[0].plot(trj2[:,0], trj2[:,1], "s")
        ax[0].plot(trj2[:,0], trj2[:,2], "s")

        ax[0].set_ylim(-0.1, 1.1)
        #ax[0].set_xlim(-0.2, 3.2)

        ax[0].plot([0,10], [0.1]*2, "--")
        ax[0].plot([0,10], [0.6]*2, "--")

        #ax[0].set_xticks(self.rep_ar[0,1:])
        #ax[1].set_xticks(self.rep_ar[0,1:])
        ax[1].set_ylim(-0.2,3.2)
        ax[1].scatter(state_df.loc[1].index.values.astype(np.float64), s_trj2.values,
                      c= state_df.loc[1]['temperature'].values, cmap = cl_cmap, s=50 )

        ax[1].plot(s_trj2, "-", lw=0.5, c='grey')
        ax[1].set_yticks([0,1,2,3])        
        return fig, ax
        
        
    def test_remd_rate_calc(self):
        state_df, all_tp_temp_df = run_multi_tp(self.remd, self.rep_ar,self.sdf)
        pt = ala4_temp_pop(state_df, self.sdf, self.rep_ar)
        rates = rates_ala4_multi_temp(all_tp_temp_df, pt, 10.0*1000.0, 4) # this factor of 10^3 should perhaps be removed from the scripts
        
        # 6 frames at 00 at temperature 1
        # 5 frames at 00 at temperature 0
        
        #during the transiton 3 frames at temperature 1, one frame at temperature 0
        r_0001_T1 = 0.75 / (10*0.6)
        r_0001_T0 = 0.25 / (10*0.5)
        r_0001 = np.array([r_0001_T0, r_0001_T1])
        
        np.testing.assert_equal( r_0001, rates.dropna()[(rates.dropna().type == (0,0,0,1))].rate)
        
        # I could test more transitions
        
    def test_find_start_end_points_tps(self):
        """
        Do I find the correct start and end points of the transition paths. Note that I include the first "stable"
        conformation in both the initial and the product state in my defintion of the transition path.
        
        See:
        /home/lustelzl/Projects/kinetics/rbn2_remd2-2/re-analysis/trj14_14jan16.png  
        """    
        #state_df, all_tp_temp_df = run_multi_tp(self.remd, self.rep_ar,self.sdf)
        _trj1, _ev1 = multi_tp(self.remd[0], self.sdf)

        events_end_d = find_end_points(_trj1, self.remd[0])
        
        np.testing.assert_equal(events_end_d[3.0], (0,0,0,1))
        np.testing.assert_equal(events_end_d[9.0], (0,1,0,0))
        
        events_start_d = find_start_points(_trj1, self.remd[0], events_end_d, self.sdf, verbose=True)
        
        np.testing.assert_equal(events_start_d[(0,0,0,1)], [(0.0, 3.0)])
        np.testing.assert_equal(events_start_d[(0,1,0,0)], [(6.0, 9.0)])
        
        
    def test_delta_T_tp(self):
        
        o = get_rates_pop_multi_temp(self.remd, self.rep_ar, "test_remd",
                                     state_def=self.sdf, return_raw_events=True)
        state_df, all_tp_temp_df, pt, rates, events_d  =o
       
        answers_d = {(1,1,0,0) : [0,1],
                     (0,0,0,1) : [0],
                     (0,1,0,0) : [0],
                     (0,0,1,1) : [0] }
        tp_trans_d =  {}
        for c in rates.type.unique():
            tp_trans_d[c] = delta_T_tps(c, events_d, self.rep_ar, trj_l=np.arange(2) )
          
        for tp, tp_delta_t in tp_trans_d.items():
            if tp_delta_t.values():
               print tp, tp_delta_t.values()
               assert answers_d[tp] == tp_delta_t.values()  
       
       
    def test_dwell_time_calculation(self):
        _s =loop_sorted_dwell_times(self.all_tp_temp_df, self.rep_ar)
        _s_gr0 = _s[_s.wait > 0]
        answers_dw_d = {}
        answers_dw_d[((0, (0,1,0,0)))] = [1.5]
        answers_dw_d[((1, (0,1,0,0)))] = [1.5]
        # traj 1 TP at 4 only involves a single TP but the wait involves two Ts
        answers_dw_d[((1, (0,0,1,1)))] = [3/ 2.0] 
        
        for key, values in answers_dw_d.items():
            assert _s_gr0[(_s_gr0.temperature == key[0] ) & (_s_gr0.type== key[1])].wait_T.values == values[0] 
            
         
    def test_tba_correction(self):
        """
        First half of a TP should belong to reactent state, second half to products.
        N.B.: If TP has odd number of frames, the additional should belong to reactents.
        How to test this?
        """

        o = get_rates_pop_multi_temp(self.mod_remd, self.rep_ar, "mod_remd",
                             state_def=self.sdf, return_raw_events=True, recalc=True)
        m_state_df, m_tp_df, m_pt, m_rates, m_events_d  =o
        m_tba = loop_corr_tba(m_state_df, m_tp_df, verbose=True)
        
        # traj 0 TP1 0,1,2,3,4
        assert np.all(m_tba.loc[0]['b'][:4].values == np.array([0,0,0,1,1]))
        print m_tba.loc[0]['b'][:4]
        
        # traj 0 TP2 6,7,8
        assert np.all(m_tba.loc[0]['b'][6:8].values == np.array([1,1,0]))
        # traj 1 TP1 0,1
        assert np.all(m_tba.loc[1]['a'][:1].values == np.array([1,0]))
        # traj 1 TP2 4,5
        assert np.all(m_tba.loc[1]['a'][4:5].values == np.array([0,1]))
        # traj 1 TP3 6,7
        assert np.all(m_tba.loc[1]['a'][6:7].values == np.array([1,0]))
        
        print  m_tba.loc[1]['a'][6:7]
        
    
    def test_av_pt_groupby_temperature(self):
        '''
        I use pandas groupby to calculate the mean populations over different simulation blocks
        See explanation:
        http://stackoverflow.com/questions/31569549/how-to-groupby-a-dataframe-in-pandas-and-keep-columns
        '''
        md1_pt = pd.read_pickle('test_data/md_st1_14apr16_bl10_pop.pickle')
        md2_pt = pd.read_pickle('test_data/md_st2_14apr16_bl10_pop.pickle')  
        md3_pt = pd.read_pickle('test_data/md_st3_14apr16_bl10_pop.pickle')
        
        _st_names = generate_state_names([[0,0],[0,0]])
        s = md1_pt[_st_names] + md2_pt[_st_names] + md3_pt[_st_names]
        av= s/ 3.0
        _pt = pd.concat([md1_pt, md2_pt, md3_pt]).groupby('temperature').mean().reset_index()
        np.testing.assert_almost_equal(_pt.values[:,1:], av.values)
        
        
    def test_combine_stages(self):
        '''
        test_combine_stages: Generate MC trajectories and compare combined rate calculation to
        manual rate calculation. 
        '''
        # MC diffusion trajectories generated on 2/5/16
        # /home/lustelzl/Projects/kinetics/prototypes/unit-test-combine-st'
        r1_st1 = np.genfromtxt('test_data/st1_mc.txt') 
        r1_st2 = np.genfromtxt('test_data/st2_mc.txt')

        mc_steps = len(r1_st1)
        t_pts = np.arange(mc_steps)
        t_pts_st2 = np.arange(mc_steps, mc_steps+mc_steps)
        
        q_st1_d = {}
        for ci, col in enumerate(r1_st1.T):
            q_st1_d[ci] = np.column_stack((t_pts, col))
            
        q_st2_d = {}
        for ci, col in enumerate(r1_st2.T):
            q_st2_d[ci] = np.column_stack((t_pts_st2, col))
        
        rep_ar_st1 = gen_rep_ar_md([0,1], sim_length=mc_steps, time_pts=t_pts)
        rep_ar_st2 = gen_rep_ar_md([0,1], sim_length=mc_steps, time_pts=t_pts_st2)
        
        rc_d = {'st1': q_st1_d, 'st2': q_st2_d}
        rep_d = {'st1' : rep_ar_st1, 'st2' : rep_ar_st2}
        
        os.chdir('test_data')
        s=[[-1, 1]]
        _c = combine_stages_df('_test',s , rc_d, rep_d,
                  stage_l=['st1', 'st2'], n_blocks=1, remd=False, tmp_l=[0,1],verbose=True,
                  time_unit_factor=1.0, dir_name='st')
        r = _c[0].dropna()
        
        # get rates from concatenated trajectories. 
        rep_comb = np.vstack((rep_ar_st1, rep_ar_st2))
        q_comb = {}
        for qi, q in q_st1_d.items():
            q_comb[qi] = np.vstack((q, q_st2_d[qi]))
        
        state_df_c, all_tp_temp_df_c, pt_c, rates_c = get_rates_pop_multi_temp(q_comb, rep_comb,
                                                                       'stacked-stages/mc_2may16',
                                                                       state_def=s , remd=False,
                                                                       time_unit_factor=1, split_state_tp=True,
                                                                       recalc=True)
                       
        r_c = rates_c.dropna()
        
        np.testing.assert_almost_equal(np.float64(r_c[r_c.type==(0,1)].rate.values),
                               np.float64(r[r.type ==(0,1)].rate.values))
                               
        np.testing.assert_almost_equal(np.float64(r_c[r_c.type==(1,0)].rate.values),
                               np.float64(r[r.type ==(1,0)].rate.values))
        
        # manual rate calculation at T=0
        # folding: 5 events
        # unfolding 6 events
        p0_u = 0.67453
        assert pt_c[(pt_c.temperature==0)].values[0,1] ==   p0_u
        
        np.testing.assert_almost_equal(5.0 / (2*mc_steps * (p0_u)),
                               r_c[(r_c.temperature == 0) & (r_c.type==(0,1))].rate)
        np.testing.assert_almost_equal(6.0 / (2*mc_steps * (1.0 -p0_u)),
                               r_c[(r_c.temperature == 0) & (r_c.type==(1,0))].rate)
                               
        # perhaps interesting test would be to run the two stages with the temperature ordering reversed- would my
        # scripts deal with this?     
        os.chdir('../') 
       
        
    def test_combine_stages_time_indeces(self):
        '''
        test_combine_stages_time_indeces: the combinations of stages should work even if the time indices
        are not set correctly in the trajectory segments (i.e. both segments to be combined start at 0)
        '''
        pass  
        
        
    def test_combine_conditions(self):
        pass

        #for tri, t in enumerate(transition_l):
        #    if tri == 0:
        #       _cond =  (tp_df.type == t)
        #    else:
        #        _cond = _cond | (tp_df.type == t)
        #pass
    
         
    def test_detect_rex_on_tp(self):
        test_stretch = [0,0,1,1,1,0,0,0,5,5,5,6]
         # 4 changes in temperature
        # lengths 2,3, 3, 3, 1
        delta_stretch_d = {}
        delta_stretch_d[(0, (0, len(test_stretch))) ] = test_stretch
        # one change in temperature
        evs = [(1, (1, 2)),  (2, (1, 1)), (3, (1, 2))]
        delta_stretch_d[evs[0]] = [0,1]
        delta_stretch_d[evs[1]] = [1]
        delta_stretch_d[evs[2]] = [4,4]
        num_rex, len_stretches = analyse_rex_on_tp(delta_stretch_d)
        
        answers_d ={(0, (0, len(test_stretch))) : [4, [2,3, 3, 3, 1]] }
        answers_d[evs[0]] = [1, [1,1]] 
        answers_d[evs[1]] = [0, [1]]
        answers_d[evs[2]] = [0, [2]]
        for ev, answer in answers_d.items():
            print ev, num_rex[ev], answer[0]
            print ev, len_stretches[ev], answer[1]
            assert num_rex[ev] == answer[0]
            assert len_stretches[ev] == answer[1]
