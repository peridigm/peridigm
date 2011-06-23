

include_directories(${UTILITIES_DIR} ${VTK_INCLUDE_DIRS})

add_subdirectory(unit_test)

# PdVTK library
add_library(PdVTK PdVTK.cxx Field.cxx)







