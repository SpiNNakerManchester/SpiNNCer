#!/bin/bash

nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.05_ms.npz \
	> spin_0.1_vs_spin_0.05.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.1_ms.npz \
	> spin_0.1_vs_spin_0.1.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.05_ms.npz results/nest_granule_test_dt_0.05_ms.npz \
	> spin_0.05_vs_nest_0.05.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/nest_granule_test_dt_0.05_ms.npz results/nest_granule_test_dt_0.1_ms.npz \
	> nest_0.05_vs_nest_0.1.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.05_ms.npz results/nest_granule_test_dt_0.1_ms.npz \
	> spin_0.05_vs_nest_0.1.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.1_ms.npz results/nest_granule_test_dt_0.1_ms.npz \
	> spin_0.1_vs_nest_0.1.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_1.0_ms.npz results/nest_granule_test_dt_1.0_ms.npz \
	> spin_1.0_vs_nest_1.0.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/nest_granule_test_dt_0.1_ms.npz results/nest_granule_test_dt_0.1_ms_ORIGINAL_SPIKES.npz \
	> nest_0.1_MODIFIED_VS_ORIGINAL_SPIKES.out 2>&1 &	

nohup python cerebellum_analysis.py --consider_delays -i results/spinn_granule_test_dt_0.05_ms.npz \
	results/spinn_granule_test_dt_0.1_ms.npz \
	results/spinn_granule_test_dt_1.0_ms.npz \
    > individual_small_spinn_tests.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays -i results/nest_granule_test_dt_0.05_ms.npz \
	results/nest_granule_test_dt_0.1_ms.npz \
	results/nest_granule_test_dt_1.0_ms.npz \
    > individual_small_nest_tests.out 2>&1 &


nohup python cerebellum_analysis.py --consider_delays -i results/nest_granule_test_dt_0.05_ms_ORIGINAL_SPIKES.npz \
	results/nest_granule_test_dt_0.1_ms_ORIGINAL_SPIKES.npz \
	results/nest_granule_test_dt_1.0_ms_ORIGINAL_SPIKES.npz \
    > individual_small_nest_tests_ORIGINAL_SPIKES.out 2>&1 &    

for i in {1..10}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.1_ms_no_loops_$i.npz results/nest_granule_test_dt_0.1_ms.npz \
	> "spin_0.1_no_loops_"$i"_vs_nest_0.1_MODIFIED_SPIKES.out" 2>&1 &
done   

for i in {1..10}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.1_ms_no_loops_$i.npz results/nest_granule_test_dt_0.1_ms_ORIGINAL_SPIKES.npz \
	> "spin_0.1_no_loops_"$i"_vs_nest_0.1_ORIGINAL_SPIKES.out" 2>&1 &
done     

# Loop at different runs 

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.1_ms_RUN_2.npz \
	> spin_0.1_vs_RUN_2.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.1_ms_no_loops_3.npz \
	> spin_0.1_vs_no_loop_3.out 2>&1 &

## CAREFUL HERE
nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file.npz \
	> full_spinn_grc_loop_3_vs_nest.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3.npz \
	> full_spinn_grc_loop_3_vs_nest_stim_3.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/spinn_from_file_400x400_3_loop.npz \
	> full_spinn_grc_loop_3_vs_ITSELF.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz \
	> full_spinn_grc_loop_3_vs_ORIGINAL_nest.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_spinn_grc_loop_3_vs_MODIFIED_SPIKES_nest_stim_3.out 2>&1 &
	

nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_ORIGINAL_SPIKES_nest_vs_MODIFIED_SPIKES_nest_stim_3.out 2>&1 &	

for i in {0..9}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  "variance_testing_3_loop_POISSON_stim_3_@0b7d8b198ae23735600f8c5e1e5cfb6f/spinn_400x400_run_"$i"_poisson_@0b7d8b198ae23735600f8c5e1e5cfb6f/results/spinn_400x400_run_"$i"_poisson_@0b7d8b198ae23735600f8c5e1e5cfb6f.npz" \
	results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> "FULL_SCALE_spin_variance_stim_3_run_"$i"_vs_nest.out" 2>&1 &
done    

# nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/nest_400x400_from_file_stim_3.npz \
# 	> full_scale_comparison_spin_0.1_vs_nest_0.1.out 2>&1 &	

# nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/spinn_400x400_run_1_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz \
# 	> full_scale_comparison_variance_testing_spin_0.1_vs_spin_0.1.out 2>&1 &	


#### SCRATCH SPACE
# nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/EXPERIMENTAL_testing_stimulus_from_file_400x400_1_loop.npz \
#  > full_spinn_vs_loop_1.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/EXPERIMENTAL_testing_stimulus_from_file_400x400_2_loop.npz \
#  > full_spinn_vs_loop_2.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz \
#  > full_spinn_vs_loop_3.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_run_0_poisson_@8adfb23bc53e1558a20d7dd2e917bf4f.npz results/EXPERIMENTAL_testing_stimulus_from_file_400x400_4_loop.npz \
#  > full_spinn_vs_loop_4.out 2>&1 &


# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_1_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_1.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_2_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_2.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_3.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_4_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_4.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_5_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_5.out 2>&1 &

# nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_6_loop.npz results/nest_400x400_from_file_stim_3.npz \
#  > full_nest_vs_loop_6.out 2>&1 &