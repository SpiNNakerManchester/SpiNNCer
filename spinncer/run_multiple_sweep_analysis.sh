# nohup python3 provenance_analysis.py -i activity_sweep_w_reinjection_f_peak_POISSON_@901769941effa128808a8ad7e3d3754c\
# 	activity_sweep_w_reinjection_f_peak_PERIODIC_@154a7bd22ba1c00611fd2fdbc8ede77f\
# 	--group_on f_peak > f_peak_analysis.out &

# nohup python3 provenance_analysis.py -i activity_sweep_w_reinjection_stim_radius_POISSON_@0a24bd8c511d32f70fe7050a37421ffb\
# 	activity_sweep_w_reinjection_stim_radius_PERIODIC_@dcbe1f3ffbd69e03195c88f6cdae5fb8\
# 	--group_on stim_radius > stim_radius_analysis.out &

# nohup python3 provenance_analysis.py -i variance_testing_from_file_@3aeb4cb341ad570b4339cafe6b0f0575 \
# 	variance_testing_@11eeb6a111f9fa1ecb947124dc70e259 \
# 	variance_testing_PERIODIC_@11eeb6a111f9fa1ecb947124dc70e259 \
# 	--group_on n_run --group_on_name run > variance_testing.out &	

# Runs with 3 loops for granule
nohup python3 provenance_analysis.py -i  \
	activity_sweep_f_peak_POISSON_@1000x_vanilla_scaling_200_grid \
	activity_sweep_f_peak_PERIODIC_@1000x_vanilla_scaling_200_grid \
	--group_on f_peak > f_peak_vanilla_grid.out &

nohup python3 provenance_analysis.py -i  \
	activity_sweep_stim_radius_POISSON_@1000x_vanilla_scaling_200_grid \
	activity_sweep_stim_radius_PERIODIC_@1000x_vanilla_scaling_200_grid \
	--group_on stim_radius > stim_radius_vanilla_grid.out &

nohup python3 provenance_analysis.py -i  
	activity_sweep_f_peak_POISSON_@1000x_vanilla_scaling_200_grid \
	activity_sweep_f_peak_3_loops_try_2_POISSON_@aea5e313e04d41712ac8f5a69031eb90 \
	activity_sweep_f_peak_3_loops_10x_PERIODIC_@547e8ae1ce1014d8ec7feed9955dccfe \
    activity_sweep_f_peak_3_loops_try_2_POISSON_@0625bb6d09db8cee6b4ca024e1e15af2 \
	activity_sweep_f_peak_3_loops_try_2_PERIODIC_@3eed4e3f66064a7f152a507b84c3f506 \
	activity_sweep_f_peak_3_loops_POISSON_@5d4e382ed4f7c3fea9f8c4d4aeece9be \
	activity_sweep_f_peak_3_loops_PERIODIC_@24c1da92e5da663c0999133e1f13132f \
	--group_on f_peak > f_peak_analysis_3_loops.out &

nohup python3 provenance_analysis.py -i activity_sweep_stim_radius_3_loops_10x_POISSON_@9e05f5aaa55d8b6711603473bf003dac \
	activity_sweep_stim_radius_3_loops_10x_PERIODIC_@46745a72771afee52c99088413c89d64 \
    activity_sweep_stim_radius_3_loops_POISSON_@9b638e178d2ac4bafa81ebd3ff9b0570 \
	activity_sweep_stim_radius_3_loops_PERIODIC_@568388bfe6a6dea2d6b35daee880837b \
	--group_on stim_radius > stim_radius_analysis_3_loops.out &

nohup python3 provenance_analysis.py -i \
	 variance_testing_3_loop_PERIODIC_@a08761636ca9e9748415c3c515695f42 \
     variance_testing_3_loop_POISSON_@0b7d8b198ae23735600f8c5e1e5cfb6f \
	 variance_testing_3_loop_POISSON_@f3503d3903819ac418e2d037cffacf85 \
	--group_on n_run --group_on_name run > variance_testing_3_loops.out &	

# Test
nohup python3 provenance_analysis.py -i \
	 variance_testing_timer_sync_no_phase_shift_POISSON_stim_3_@0e5b02b0603e04165747acce37aecc83 \
	--group_on n_run --group_on_name run > variance_testing_3_loops_no_phase_shift_timer_sync.out &	

	# Same board test
nohup python3 ../provenance_analysis.py -i \
	 variance_testing_POISSON_stim_3_@9ffa4b1557f97904a92d0db09a3c25fb \
	--group_on n_run --group_on_name run > variance_testing_same_board.out &	

# Same board test
nohup python3 provenance_analysis.py -i \
	 variance_testing_timer_sync_with_phase_shift_POISSON_stim_3_@771ddf9a9f2dd54d609a57cb4e65ea22 \
	--group_on n_run --group_on_name run > variance_testing_timer_sync.out &	
	