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

nohup python cerebellum_analysis.py --consider_delays --compare  results/nest_granule_test_dt_0.05_ms.npz results/nest_granule_test_dt_0.05_ms_ORIGINAL_SPIKES.npz \
	> nest_0.05_MODIFIED_VS_ORIGINAL_SPIKES.out 2>&1 &	

nohup python cerebellum_analysis.py --consider_delays --compare  results/nest_granule_test_dt_1.0_ms.npz results/nest_granule_test_dt_1.0_ms_ORIGINAL_SPIKES.npz \
	> nest_1.0_MODIFIED_VS_ORIGINAL_SPIKES.out 2>&1 &	

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
	nohup python cerebellum_analysis.py --consider_delays -i  results/spinn_granule_test_dt_0.1_ms_no_loops_$i.npz \
	> "INDIVIDUAL_grc_test_spin_0.1_no_loops_"$i"_.out" 2>&1 &
done  

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

nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_run_7_poisson_@9770e92909db6157c658a63a1c93c33d.npz results/nest_400x400_from_file_stim_3.npz \
	> full_spinn_grc_loop_3_vs_nest_stim_3.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/spinn_from_file_400x400_3_loop.npz \
	> full_spinn_grc_loop_3_vs_ITSELF.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz \
	> full_spinn_grc_loop_3_vs_ORIGINAL_nest.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/EXPERIMENTAL_testing_stimulus_from_file_400x400_3_loop.npz results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_spinn_grc_loop_3_vs_MODIFIED_SPIKES_nest_stim_3.out 2>&1 &
	

nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_ORIGINAL_SPIKES_nest_vs_MODIFIED_SPIKES_nest_stim_3.out 2>&1 &	

nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_from_file_stim_3.npz results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_nest_vs_MODIFIED_SPIKES_nest_stim_3.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz results/nest_400x400_from_file_stim_3.npz \
	> full_ORIGINAL_SPIKES_nest_v_nest_stim_3.out 2>&1 &	

for i in {0..9}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  "variance_testing_3_loop_POISSON_stim_3_@9770e92909db6157c658a63a1c93c33d/spinn_400x400_run_"$i"_poisson_@9770e92909db6157c658a63a1c93c33d/results/spinn_400x400_run_"$i"_poisson_@9770e92909db6157c658a63a1c93c33d.npz" \
	results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> "FULL_SCALE_spin_variance_stim_3_run_"$i"_vs_nest_MODIFIED_SPIKES.out" 2>&1 &
done    

for i in {0..9}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  "variance_testing_3_loop_POISSON_stim_3_@9770e92909db6157c658a63a1c93c33d/spinn_400x400_run_"$i"_poisson_@9770e92909db6157c658a63a1c93c33d/results/spinn_400x400_run_"$i"_poisson_@9770e92909db6157c658a63a1c93c33d.npz" \
	results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz \
	> "FULL_SCALE_spin_variance_stim_3_run_"$i"_vs_nest_ORIGINAL_SPIKES.out" 2>&1 &
done  


# TEST RUNS WITH TIMER SYNC -- WINNER
for i in {0..5}
do 
	let j=$i+1
	file_1="spinn_400x400_run_${i}_poisson_@9ffa4b1557f97904a92d0db09a3c25fb"
	file_2="spinn_400x400_run_${j}_poisson_@9ffa4b1557f97904a92d0db09a3c25fb"
	nohup python ../cerebellum_analysis.py --consider_delays --compare \
	"variance_testing_POISSON_stim_3_@9ffa4b1557f97904a92d0db09a3c25fb/"$file_1"/results/"$file_1".npz" \
	"variance_testing_POISSON_stim_3_@9ffa4b1557f97904a92d0db09a3c25fb/"$file_2"/results/"$file_2".npz" \
	> "FULL_SCALE_timer_sync_spin_variance_stim_3_run_"$i"_vs_run_"$j".out" 2>&1 &
done  

# TEST RUNS WITH TIMER SYNC -- SAME BOARD
for i in {0..8}
do 
	let j=$i+1
	file_1="spinn_400x400_run_${i}_poisson_@771ddf9a9f2dd54d609a57cb4e65ea22"
	file_2="spinn_400x400_run_${j}_poisson_@771ddf9a9f2dd54d609a57cb4e65ea22"
	nohup python cerebellum_analysis.py --consider_delays --compare \
	"variance_testing_timer_sync_with_phase_shift_POISSON_stim_3_@771ddf9a9f2dd54d609a57cb4e65ea22/"$file_1"/results/"$file_1".npz" \
	"variance_testing_timer_sync_with_phase_shift_POISSON_stim_3_@771ddf9a9f2dd54d609a57cb4e65ea22/"$file_2"/results/"$file_2".npz" \
	> "FULL_SCALE_timer_sync_spin_variance_stim_3_run_"$i"_vs_run_"$j".out" 2>&1 &
done  

# TEST RUNS WITH TIMER SYNC VS NEST
for i in {0..8}
do 
	let j=$i+1
	file_1="spinn_400x400_run_${i}_poisson_@771ddf9a9f2dd54d609a57cb4e65ea22"
	file_2="spinn_400x400_run_${j}_poisson_@771ddf9a9f2dd54d609a57cb4e65ea22"
	nohup python cerebellum_analysis.py --consider_delays --compare \
	"variance_testing_timer_sync_with_phase_shift_POISSON_stim_3_@771ddf9a9f2dd54d609a57cb4e65ea22/"$file_1"/results/"$file_1".npz" \
	results/nest_400x400_from_file_stim_3.npz \
	> "FULL_SCALE_timer_sync_spin_variance_stim_3_run_"$i"_vs_NEST.out" 2>&1 &
done 

# TEST RUNS WITHOUT TIMER SYNC -- LOOSER
for i in {0..8}
do 
	let j=$i+1
	suff_1="1e16a23b66a669dcfc9f2bfe1759eab4"
	suff_2="1e16a23b66a669dcfc9f2bfe1759eab4"
	file_1="spinn_400x400_run_${i}_poisson_@${suff_1}"
	file_2="spinn_400x400_run_${j}_poisson_@${suff_2}"
	nohup python ../cerebellum_analysis.py --consider_delays --compare \
	"variance_testing_POISSON_stim_3_@${suff_1}/"$file_1"/results/"$file_1".npz" \
	"variance_testing_POISSON_stim_3_@${suff_2}/"$file_2"/results/"$file_2".npz" \
	> "FULL_SCALE_timer_sync_spin_variance_stim_3_run_"$i"_vs_run_"$j".out" 2>&1 &
done  


nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_from_file_stim_3.npz \
	> full_nest_stim_3.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_from_file_stim_3_MODIFIED_SPIKES.npz \
	> full_nest_MODIFIED_SPIKES_stim_3.out 2>&1 &	
nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_from_file_stim_3_ORIGINAL_SPIKES.npz \
	> full_nest_ORIGINAL_SPIKES_stim_3.out 2>&1 &	   


nohup python cerebellum_analysis.py --consider_delays --i results/spinn_400x400_run_7_poisson_@9770e92909db6157c658a63a1c93c33d.npz  \
	> INDIVIDUAL_full_3_loops_analysis.out 2>&1 &

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