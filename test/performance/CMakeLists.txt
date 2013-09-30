#
# Add the tests
#
if(TEST_PERFORMANCE)
  add_test (disk_impact_np4 python ./disk_impact/disk_impact.py -machine ${MACHINE})
endif(TEST_PERFORMANCE)

add_custom_target( ptest
   COMMAND ctest
)
