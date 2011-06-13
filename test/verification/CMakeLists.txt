# List tests to be run by "make test"
#
# Create symbolic links to files in Test directory
#
#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for verification tests for out-of-source build\n")
   #
   # Loop over all files in the directory tree
   # Copy directory structure and create symbolic links to files
   #
#   file (GLOB_RECURSE directoryContents "${CMAKE_SOURCE_DIR}/test/verification/*")
#   foreach (item ${directoryContents})
      #
      # Create the subdirectory, if needed
      #
#      get_filename_component (itemPath ${item} PATH)
#      file (RELATIVE_PATH relPath "${CMAKE_SOURCE_DIR}/test/verification" ${itemPath})
#      file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test/verification/${relPath})
	  #
      # Create symbolic link to file
      #
#      get_filename_component (itemName ${item} NAME)
#      execute_process (COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/test/verification/${relPath}/${itemName} ${CMAKE_BINARY_DIR}/test/verification/${relPath}/${itemName})  
#   endforeach (item)
#endif (NOT ${buildInSource})
#
# Add the tests
#
add_test (Compression_2x1x1_np1 python ./Compression_2x1x1/np1/Compression_2x1x1.py)
add_test (Compression_2x1x1_np2 python ./Compression_2x1x1/np2/Compression_2x1x1.py)
add_test (CompressionInitDisp_2x1x1_np1 python ./CompressionInitDisp_2x1x1/np1/CompressionInitDisp_2x1x1.py)
add_test (CompressionInitDisp_2x1x1_np2 python ./CompressionInitDisp_2x1x1/np2/CompressionInitDisp_2x1x1.py)
add_test (CompressionImplicit_2x1x1_np1 python ./CompressionImplicit_2x1x1/np1/CompressionImplicit_2x1x1.py)
add_test (CompressionImplicitEssentialBC_2x1x1_np1 python ./CompressionImplicitEssentialBC_2x1x1/np1/CompressionImplicitEssentialBC_2x1x1.py)
add_test (Compression_3x1x1_np1 python ./Compression_3x1x1/np1/Compression_3x1x1.py)
add_test (Compression_3x1x1_np2 python ./Compression_3x1x1/np2/Compression_3x1x1.py)
add_test (Compression_3x1x1_np3 python ./Compression_3x1x1/np3/Compression_3x1x1.py)
add_test (BondBreakingInitialVelocity_np1 python ./BondBreakingInitialVelocity/np1/BondBreakingInitialVelocity.py)
add_test (BondBreakingInitialVelocity_np3 python ./BondBreakingInitialVelocity/np3/BondBreakingInitialVelocity.py)
add_test (BondBreakingInitialVelocity-EP_np1 python ./BondBreakingInitialVelocity-EP/np1/BondBreakingInitialVelocity.py)
add_test (BondBreakingInitialVelocity-EP_np3 python ./BondBreakingInitialVelocity-EP/np3/BondBreakingInitialVelocity.py)
add_test (Contact_2x1x1_np1 python ./Contact_2x1x1/np1/Contact_2x1x1.py)
add_test (Contact_2x1x1_np2 python ./Contact_2x1x1/np2/Contact_2x1x1.py)
if (USE_DAKOTA)
   add_test (Dakota_2x1x1_np1 python ./Dakota_2x1x1/np1/Dakota_2x1x1.py)
endif(USE_DAKOTA)
#add_test (Contact_2x1x1_np1 python ./Contact_2x1x1/np1/Contact_2x1x1.py)

add_custom_target( vtest
   COMMAND ctest
)
