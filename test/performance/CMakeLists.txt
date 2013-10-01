#
# Add the tests
#
if(PERFORMANCE_TEST_MACHINE)
  add_test (disk_impact_performance_np4 python ./disk_impact/disk_impact.py -machine ${PERFORMANCE_TEST_MACHINE})
endif()

add_custom_target( ptest
   COMMAND ctest
)
