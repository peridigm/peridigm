add_subdirectory (unitTest/)

#
# Custom Commands
#
file (GLOB Peridigm_COMPUTE_HEADERS "Peridigm*.hpp")

add_custom_command(
   OUTPUT compute_includes.hpp
   DEPENDS ${Peridigm_COMPUTE_HEADERS}
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/compute
   COMMAND sh ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.sh
)
ADD_CUSTOM_TARGET(generate_compute_includes DEPENDS compute_includes.hpp)
