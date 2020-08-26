#!/bin/bash

nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.05_ms.npz \
	> spin_0.1_vs_spin_0.05.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_granule_test_dt_0.1_ms.npz results/spinn_granule_test_dt_0.1_ms.npz \
	> spin_0.1_vs_spin_0.1.out 2>&1 &

nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_granule_test_dt_0.05_ms.npz results/nest_granule_test_dt_0.05_ms.npz \
	> spin_0.05_vs_nest_0.05.out 2>&1 &

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

nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_pss_3_rtn_accum_vanilla.npz results/nest_400x400_pss_3_rtn_accum_r_mem.npz \
	> nest_vanilla_vs_rmem.out 2>&1 &
nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_pss_3_rtn_accum_vanilla.npz results/nest_400x400_pss_3.npz \
> nest_vanilla_vs_original.out 2>&1 &	


# Another look at nest values

nohup python cerebellum_analysis.py --worst_case_spikes --consider_delays -i results/nest_400x400_pss_3.npz \
	> full_nest_stim_3.out 2>&1 &	


nohup python cerebellum_analysis.py --worst_case_spikes --consider_delays -i results/nest_400x400_pss_seed.npz \
	> full_nest_seed.out 2>&1 &		


nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_retest_400x400_poisson_r_mem_stim_3_stochastic_rounding.npz  results/nest_400x400_pss_3.npz \
	> stochastic_round_vs_nest_stim_3.out 2>&1 &	


nohup python cerebellum_analysis.py --consider_delays --compare results/spinn_retest_400x400_poisson_vanilla_stim_3_stochastic_rounding.npz  results/nest_400x400_pss_3.npz \
	> vanilla_stochastic_round_vs_nest_stim_3.out 2>&1 &

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


# R_mem EXPERIMENTS

for i in {0..9}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  \
	"test_same_board/variance_testing_POISSON_@1000x_r_mem_test_stim_3/spinn_400x400_run_"$i"_poisson_@1000x_r_mem_test_stim_3/results/spinn_400x400_run_"$i"_poisson_@1000x_r_mem_test_stim_3.npz" \
	results/nest_400x400_pss_3.npz \
	> "1000x_r_mem_spin_variance_stim_3_run_"$i"_vs_nest.out" 2>&1 &
done  

for i in {0..9}
do 
	nohup python cerebellum_analysis.py --consider_delays --compare  \
	"test_same_board_100x/variance_testing_POISSON_@100x_r_mem_test_stim_3/spinn_400x400_run_"$i"_poisson_@100x_r_mem_test_stim_3/results/spinn_400x400_run_"$i"_poisson_@100x_r_mem_test_stim_3.npz" \
	results/nest_400x400_pss_3.npz \
	> "100x_r_mem_spin_variance_stim_3_run_"$i"_vs_nest.out" 2>&1 &
done  

for i in {0..4}
do 
	nohup python ../../cerebellum_analysis.py --consider_delays --compare  \
	"spinn_400x400_run_"$i"_poisson_@1000x_scaling_200_grid_SR_r_mem_loops_4/results/spinn_400x400_run_"$i"_poisson_@1000x_scaling_200_grid_SR_r_mem_loops_4.npz" \
	../../results/nest_400x400_pss_3.npz \
	> "100x_r_mem_spin_variance_stim_3_run_"$i"_vs_nest.out" 2>&1 &
done  

# Test results from a single board with R_mem
for i in {0..4}
do 
	let j=$i+1
	file_1="spinn_400x400_run_${i}_poisson_@1000x_r_mem_scale_200_stim_3"
	file_2="spinn_400x400_run_${j}_poisson_@1000x_r_mem_scale_200_stim_3"
	nohup python ../../cerebellum_analysis.py --consider_delays --compare \
	$file_1"/results/"$file_1".npz" \
	$file_2"/results/"$file_2".npz" \
	> "1000x_r_mem_test_stim_3_run_"$i"_vs_run_"$j".out" 2>&1 &
done  

for i in {0..8}
do 
	let j=$i+1
	file_1="spinn_400x400_run_${i}_poisson_@100x_r_mem_test_stim_3"
	file_2="spinn_400x400_run_${j}_poisson_@100x_r_mem_test_stim_3"
	nohup python cerebellum_analysis.py --consider_delays --compare \
	"test_same_board_100x/variance_testing_POISSON_@100x_r_mem_test_stim_3/"$file_1"/results/"$file_1".npz" \
	"test_same_board_100x/variance_testing_POISSON_@100x_r_mem_test_stim_3/"$file_2"/results/"$file_2".npz" \
	> "100x_r_mem_test_stim_3_run_"$i"_vs_run_"$j".out" 2>&1 &
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
	file_1="spinn_400x400_run_${i}_poisson_@test_determinism_spinn_tools"
	file_2="spinn_400x400_run_${j}_poisson_@test_determinism_spinn_tools"
	nohup python cerebellum_analysis.py --consider_delays --compare \
	"test_same_board_100x/variance_testing_POISSON_@test_determinism_spinn_tools/"$file_1"/results/"$file_1".npz" \
	"test_same_board_100x/variance_testing_POISSON_@test_determinism_spinn_tools/"$file_2"/results/"$file_2".npz" \
	> "100x_comparison_queue_user_event_run_"$i"_vs_run_"$j".out" 2>&1 &
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

for i in {1..5}
do 
nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_pss_3_SR_vanilla_no_loops_$i.npz results/nest_400x400_pss_3.npz \
> "spinn_400x400_pss_3_SR_VANILLA_no_loops_"$i"_vs_nest.out" 2>&1 & 
done

for i in {1..5}
do 
nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_pss_3_SR_vanilla_no_loops_$i.npz results/nest_400x400_pss_3_rtn_accum_vanilla.npz \
> "spinn_400x400_pss_3_SR_VANILLA_no_loops_"$i"_vs_nest_rtn_accum_VANILLA.out" 2>&1 & 
done

for i in {1..5}
do 
nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_pss_3_SR_vanilla_no_loops_$i.npz results/nest_400x400_pss_3_rtn_accum_vanilla.npz \
> "spinn_400x400_pss_3_SR_VANILLA_no_loops_"$i"_vs_nest_rtn_accum_VANILLA.out" 2>&1 & 
done

for i in {1..5}
do 
nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_pss_3_SR_vanilla_no_loops_$i.npz results/spinn_400x400_pss_3_SR_r_mem_no_loops_$i.npz \
> "spinn_400x400_pss_3_SR_VANILLA_vs_R_mem_no_loops_"$i".out" 2>&1 & 
done

# if I set NEST up so that it uses sPyNNaker weights 
for i in {1..5}
do 
nohup python cerebellum_analysis.py --consider_delays --compare  results/spinn_400x400_pss_3_SR_vanilla_no_loops_$i.npz \
results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
> "spinn_400x400_pss_3_SR_VANILLA_no_loops_"$i"_vs_nest_rtn_accum_weights_and_IO_VANILLA.out" 2>&1 & 
done



#### Extra tests with RTN i_offset
nohup python3 cerebellum_experiment.py  -o spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset \
-i scaffold_full_dcn_400.0x400.0_v3.hdf5 -s 400x400_stimulus_3.npz > spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset.out 2>&1 &

nohup python3 ../cerebellum_experiment.py  -o spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset \
-i scaffold_full_dcn_400.0x400.0_v3.hdf5 -s 400x400_stimulus_3.npz > spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset.out 2>&1 &


# SR vs NEST
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
    > "spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset_vs_nest_IO_RTN.out" 2>&1 &

nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
    > "spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset_vs_nest_IO_RTN.out" 2>&1 &

nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss.npz \
    > "spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset_vs_nest.out" 2>&1 &

nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss.npz \
    > "spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset_vs_nest.out" 2>&1 &    

# RTN vs NEST
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_RTN_VANILLA_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
    > "spinn_400x400_RTN_VANILLA_pss_3_rtn_ioffset_vs_nest_IO_RTN.out" 2>&1 &

nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_RTN_R_MEM_pss_3_rtn_ioffset.npz  ../results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
    > "spinn_400x400_RTN_R_MEM_pss_3_rtn_ioffset_vs_nest_IO_RTN.out" 2>&1 &
# SR vs SR
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset.npz  results/spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset.npz \
    > "spinn_400x400_SR_VANILLA_vs_R_MEM.out" 2>&1 &
# SR vs RTN
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_VANILLA_pss_3_rtn_ioffset.npz  results/spinn_400x400_RTN_VANILLA_pss_3_rtn_ioffset.npz \
    > "spinn_400x400_VANILLA_SR_vs_RTN.out" 2>&1 &
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_SR_R_MEM_pss_3_rtn_ioffset.npz  results/spinn_400x400_RTN_R_MEM_pss_3_rtn_ioffset.npz \
    > "spinn_400x400_R_MEM_SR_vs_RTN.out" 2>&1 &

# RTN vs RTN
nohup python ../cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_RTN_VANILLA_pss_3_rtn_ioffset.npz  results/spinn_400x400_RTN_R_MEM_pss_3_rtn_ioffset.npz \
    > "spinn_400x400_RTN_VANILLA_vs_R_MEM.out" 2>&1 &

#### NEST AREA 


nohup python cerebellum_analysis.py --consider_delays --compare results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz  results/nest_400x400_pss_3.npz \
    > "nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset_vs_full_accuracy.out" 2>&1 & 

nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.npz \
    > "an_nest_400x400_pss_3_rtn_accum_vanilla_accum_ioffset.out" 2>&1 & 

nohup python cerebellum_analysis.py --consider_delays --worst_case_spikes -i results/nest_400x400_periodic_retest.npz \
    > "an_nest_400x400_periodic_retest.out" 2>&1 & 

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