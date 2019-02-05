if(HAVE_YAML)
  add_test (Notched_Plate_Data_Loader_np1 python ./np1/Notched_Plate_Data_Loader.py)
  add_test (Notched_Plate_Data_Loader_np4 python ./np4/Notched_Plate_Data_Loader.py)
endif(HAVE_YAML)
