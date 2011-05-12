file (GLOB computeSOURCES 
  "${CMAKE_SOURCE_DIR}/src/compute/*.cpp"
)

# Tell CMake a header file won't be available until build time.
set_source_files_properties(./compute/compute_includes.hpp PROPERTIES GENERATED 1)

set(utPeridigm_Compute_Force_SOURCES
    ../../Peridigm_DiscretizationFactory.cpp
    ../../Peridigm_PdQuickGridDiscretization.cpp
    ../../Peridigm_DataManager.cpp
    ../../Peridigm_State.cpp
    ../../Peridigm_ComputeManager.cpp
    ${computeSOURCES}
    ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
)
if(PERIDIGM_STK)
  set(utPeridigm_Compute_Force_SOURCES ../../Peridigm_STKDiscretization.cpp ${utPeridigm_Compute_Force_SOURCES})
endif()

add_executable(utPeridigm_Compute_Force ./utPeridigm_Compute_Force.cpp ${utPeridigm_Compute_Force_SOURCES}) 
target_link_libraries(utPeridigm_Compute_Force ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES} )

#
# Custom Commands
#
file (GLOB Peridigm_COMPUTE_HEADERS
"${CMAKE_SOURCE_DIR}/src/compute/Peridigm*.hpp")

add_custom_command(
   OUTPUT  ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
   DEPENDS ${Peridigm_COMPUTE_HEADERS}
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/compute
   COMMAND sh ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.sh
)
