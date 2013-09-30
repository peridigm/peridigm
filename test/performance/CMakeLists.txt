#
# Add the tests
#
if(TEST_PERFORMANCE)
  add_test (elastic_disk_impact_np4 python ./elastic_disk_impact/elastic_disk_impact.py)
endif(TEST_PERFORMANCE)

add_custom_target( ptest
   COMMAND ctest
)
