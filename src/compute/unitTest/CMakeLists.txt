file (GLOB computeSOURCES "${CMAKE_SOURCE_DIR}/src/compute/*.cpp")

set(utPeridigm_Compute_Force_SOURCES
    ../../Peridigm_DiscretizationFactory.cpp
    ../../Peridigm_PdQuickGridDiscretization.cpp
    ../../Peridigm_DataManager.cpp
    ../../Peridigm_State.cpp
    ../../Peridigm_ComputeManager.cpp
    ${computeSOURCES}
)
if(PERIDIGM_STK)
  set(utPeridigm_Compute_Force_SOURCES ../../Peridigm_STKDiscretization.cpp ${utPeridigm_Compute_Force_SOURCES})
endif()

add_executable(utPeridigm_Compute_Force ./utPeridigm_Compute_Force.cpp ${utPeridigm_Compute_Force_SOURCES}) 
target_link_libraries(utPeridigm_Compute_Force ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES} )
add_dependencies(utPeridigm_Compute_Force generate_compute_includes)
