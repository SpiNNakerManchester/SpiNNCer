nohup python cerebellum_analysis.py --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_140_matching.npz > analysis_gold_400_stim_radius_140_matching.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_70.npz > analysis_gold_400_stim_radius_70.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/gold_standards/gold_standard_results_158.npz > analysis_gold_158.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_140.npz > analysis_gold_400_stim_radius_140.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/gold_standards/gold_standard_results_400_fcn.npz > analysis_gold_standard_results_400_fcn.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/nest_400x400.npz > analysis_nest_400x400.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/nest_158x158.npz > analysis_nest_158x158.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i 140um_stim_experiments_2000x_pss/results/2000x_pss_400x400_140um_re.npz > analysis_2000x_pss_400x400_140um_re.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i 140um_stim_experiments_2000x/results/2000x_ssa_400x400_140um_re.npz > analysis_2000x_ssa_400x400_140um_re.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/nest_400x400_no_proj.npz > analysis_nest_400x400_no_proj.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/nest_158x158_no_proj.npz > analysis_nest_158x158_no_proj.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/spinnaker_158x158.npz > analysis_spinnaker_158x158.out 2>&1 &
nohup python cerebellum_analysis.py --worst_case_spikes -i results/spinnaker_158x158_no_proj.npz > analysis_spinnaker_158x158_no_proj.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_140_matching.npz > analysis_gold_400_stim_radius_140_matching_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_70.npz > analysis_gold_400_stim_radius_70_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/gold_standards/gold_standard_results_158.npz > analysis_gold_158_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/gold_standards/gold_standard_results_400_stim_radius_140.npz > analysis_gold_400_stim_radius_140_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/gold_standards/gold_standard_results_400_fcn.npz > analysis_gold_standard_results_400_fcn_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400.npz > analysis_nest_400x400_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_158x158.npz > analysis_nest_158x158_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i 140um_stim_experiments_2000x_pss/results/2000x_pss_400x400_140um_re.npz > analysis_2000x_pss_400x400_140um_re_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i 140um_stim_experiments_2000x/results/2000x_ssa_400x400_140um_re.npz > analysis_2000x_ssa_400x400_140um_re_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_no_proj.npz > analysis_nest_400x400_no_proj_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_158x158_no_proj.npz > analysis_nest_158x158_no_proj_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/spinnaker_158x158.npz > analysis_spinnaker_158x158_w_delays.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/spinnaker_158x158_no_proj.npz > analysis_spinnaker_158x158_no_proj.out 2>&1 &
