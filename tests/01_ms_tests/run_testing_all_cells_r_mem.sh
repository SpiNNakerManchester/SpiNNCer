rounding=RTN
nohup python testing_all_cells_single_spikes.py --simtime 10000 --suffix "$rounding" --r_mem > spinn_R_MEM_"$rounding"_single_spike.out  2>&1 &
nohup python testing_all_cells_spikes_realistic_conn.py --simtime 10000  --suffix "$rounding" --r_mem > spinn_R_MEM_"$rounding"_real_conn_test.out  2>&1 &
nohup python testing_all_cells_dc.py --simtime 10000  --suffix "$rounding" --r_mem > spinn_R_MEM_"$rounding"_DC.out  2>&1 &