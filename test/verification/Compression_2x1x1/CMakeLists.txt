if (HAVE_YAML)
   add_test (Compression_2x1x1_np1 python ./np1/Compression_2x1x1.py)
   add_test (Compression_2x1x1_np2 python ./np2/Compression_2x1x1.py)
endif (HAVE_YAML)
