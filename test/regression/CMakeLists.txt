# List tests to be run by "make test"
add_test (RebalanceCube python ./RebalanceCube/RebalanceCube.py)

add_custom_target( rtest
   COMMAND ctest
)
