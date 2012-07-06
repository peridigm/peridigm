#
# Add the tests
#
#add_test (Compression_2x1x1_np1 python ./Compression_2x1x1/np1/Compression_2x1x1.py)
#add_test (Compression_2x1x1_np2 python ./Compression_2x1x1/np2/Compression_2x1x1.py)
#add_test (CompressionInitDisp_2x1x1_np1 python ./CompressionInitDisp_2x1x1/np1/CompressionInitDisp_2x1x1.py)
#add_test (CompressionInitDisp_2x1x1_np2 python ./CompressionInitDisp_2x1x1/np2/CompressionInitDisp_2x1x1.py)
#add_test (CompressionImplicit_2x1x1_np1 python ./CompressionImplicit_2x1x1/np1/CompressionImplicit_2x1x1.py)
#####add_test (CompressionImplicit_2x1x1_np2 python ./CompressionImplicit_2x1x1/np2/CompressionImplicit_2x1x1.py)
#add_test (CompressionImplicitEssentialBC_2x1x1_np1 python ./CompressionImplicitEssentialBC_2x1x1/np1/CompressionImplicitEssentialBC_2x1x1.py)
#####add_test (CompressionImplicitEssentialBC_2x1x1_np2 python ./CompressionImplicitEssentialBC_2x1x1/np2/CompressionImplicitEssentialBC_2x1x1.py)
#add_test (Compression_3x1x1_np1 python ./Compression_3x1x1/np1/Compression_3x1x1.py)
#add_test (Compression_3x1x1_np2 python ./Compression_3x1x1/np2/Compression_3x1x1.py)
#add_test (Compression_3x1x1_np3 python ./Compression_3x1x1/np3/Compression_3x1x1.py)
#add_test (BondBreakingInitialVelocity_np1 python ./BondBreakingInitialVelocity/np1/BondBreakingInitialVelocity.py)
#add_test (BondBreakingInitialVelocity_np3 python ./BondBreakingInitialVelocity/np3/BondBreakingInitialVelocity.py)
#add_test (BondBreakingInitialVelocity-EP_np1 python ./BondBreakingInitialVelocity-EP/np1/BondBreakingInitialVelocity.py)
#add_test (BondBreakingInitialVelocity-EP_np3 python ./BondBreakingInitialVelocity-EP/np3/BondBreakingInitialVelocity.py)
add_test (Contact_2x1x1_np1 python ./Contact_2x1x1/np1/Contact_2x1x1.py)
add_test (Contact_2x1x1_np2 python ./Contact_2x1x1/np2/Contact_2x1x1.py)
#if (USE_DAKOTA)
#   add_test (Dakota_2x1x1_np1 python ./Dakota_2x1x1/np1/Dakota_2x1x1.py)
#endif(USE_DAKOTA)
#####add_test (ep_cube_np1 python ./ep_cube/np1/ep_cube_np1.py)
#####add_test (ep_cube_np2 python ./ep_cube/np2/ep_cube_np2.py)

add_custom_target( vtest
   COMMAND ctest
)
