#
# Add the tests
#
if(PERFORMANCE_TEST_MACHINE)
  add_test (disk_impact_performance_np4 python ./disk_impact/disk_impact.py -machine ${PERFORMANCE_TEST_MACHINE})
  add_test (tensile_test_performance_np4 python ./tensile_test/tensile_test.py -machine ${PERFORMANCE_TEST_MACHINE})
  add_test (twist_and_pull_performance_np4 python ./twist_and_pull/twist_and_pull.py -machine ${PERFORMANCE_TEST_MACHINE})
  add_test (fragmenting_cylinder_performance_np4 python ./fragmenting_cylinder/fragmenting_cylinder.py -machine ${PERFORMANCE_TEST_MACHINE})
  #
  # this "test" is actually a script that parses the performance test results and writes a table to test/performance/test_summary/performance_test_summary.txt
  #
  add_test (summarize_performace_tests python ./test_summary/query_performance_test_results.py)
endif()

add_custom_target( ptest
   COMMAND ctest
)
