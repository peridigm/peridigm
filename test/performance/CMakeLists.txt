#
# Add the tests
#
if(PERFORMANCE_TEST_MACHINE)
  add_test (disk_impact_performance_np4 python ./disk_impact/disk_impact.py -machine ${PERFORMANCE_TEST_MACHINE})
  add_test (tensile_test_performance_np4 python ./tensile_test/tensile_test.py -machine ${PERFORMANCE_TEST_MACHINE})
  add_test (twist_and_pull_performance_np4 python ./twist_and_pull/twist_and_pull.py -machine ${PERFORMANCE_TEST_MACHINE})
endif()

add_custom_target( ptest
   COMMAND ctest
)
