include_directories(${UTILITIES_DIR} ${VTK_INCLUDE_DIRS})

# PdField library
add_library(PdField Field.cxx)
