## Makefile.am -- Process this file with automake to produce Makefile.in

######################################################################
# EXECUTOR
######################################################################

check_PROGRAMS += \
				  logical_tile_test \
				  mutate_test \
				  materialization_test \
				  seq_scan_test \
				  index_scan_test \
				  limit_test \
				  join_test \
				  order_by_test \
				  hash_set_op_test \
				  aggregate_test \
				  append_test \
				  projection_test \
				  tile_group_layout_test

executor_tests_common= 	executor/executor_tests_util.cpp \
						harness.cpp

logical_tile_test_SOURCES = \
							$(executor_tests_common) \
							executor/logical_tile_test.cpp 

materialization_test_SOURCES = \
							   $(executor_tests_common) \
							   executor/materialization_test.cpp 

mutate_test_SOURCES = \
					  $(executor_tests_common) \
					  executor/mutate_test.cpp 

seq_scan_test_SOURCES = \
						$(executor_tests_common) \
						executor/seq_scan_test.cpp 

index_scan_test_SOURCES = \
						  $(executor_tests_common) \
						  executor/index_scan_test.cpp 

limit_test_SOURCES = \
					 $(executor_tests_common) \
					 executor/limit_test.cpp
								
join_test_SOURCES = \
						$(executor_tests_common) \
						executor/join_test.cpp \
						executor/join_tests_util.cpp

order_by_test_SOURCES = \
						$(executor_tests_common) \
						executor/order_by_test.cpp 
					
hash_set_op_test_SOURCES = \
						$(executor_tests_common) \
						executor/hash_set_op_test.cpp
						
aggregate_test_SOURCES = \
						$(executor_tests_common) \
						executor/aggregate_test.cpp 
					
append_test_SOURCES = \
					$(executor_tests_common) \
					executor/append_test.cpp					

projection_test_SOURCES = \
					$(executor_tests_common) \
					executor/projection_test.cpp					
					
tile_group_layout_test_SOURCES = \
								 $(executor_tests_common) \
								 executor/tile_group_layout_test.cpp
