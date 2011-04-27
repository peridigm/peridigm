# List tests to be run by "make test"
add_test (RebalanceCube python ./RebalanceCube/RebalanceCube.py)
add_test (Compression_QS_3x2x2_np1 python ./Compression_QS_3x2x2/np1/Compression_QS_3x2x2.py)
add_test (Compression_QS_3x2x2_np2 python ./Compression_QS_3x2x2/np2/Compression_QS_3x2x2.py)

add_custom_target( rtest
   COMMAND ctest
)
