# global_placement -use_smt_threads runs every requested thread instead of
# capping placement at the physical core count (GPL-0190). The log has no
# GPL-0190 line and the placement matches core01's single-threaded golden.
source helpers.tcl
set test_name smt_threads01
read_lef ./nangate45.lef
read_def ./core01.def

set_thread_count 16
global_placement -density 0.6 -init_density_penalty 0.01 -skip_initial_place \
  -use_smt_threads

set def_file [make_result_file $test_name.def]
write_def $def_file
diff_files $def_file core01.defok
