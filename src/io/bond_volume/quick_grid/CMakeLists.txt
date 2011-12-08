#
# QuickGrid library
#
add_subdirectory(unit_test)
include_directories(${CMAKE_SOURCE_DIR}/io ${CMAKE_SOURCE_DIR})
add_library(QuickGridVolFrac calculators.cxx)
