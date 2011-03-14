#
# QuickGrid Unit Tests
#

include_directories(${UTILITIES_DIR})

add_executable(ut_Q2CylinderRingHorizon ut_Q2CylinderRingHorizon.cxx)
target_link_libraries(ut_Q2CylinderRingHorizon  QuickGrid ${UT_REQUIRED_LIBS})
add_test (ut_Q2CylinderRingHorizon python run_unit_test.py ./ut_Q2CylinderRingHorizon)

add_executable(ut_QuickGridHorizon ut_QuickGridHorizon.cxx)
target_link_libraries(ut_QuickGridHorizon  QuickGrid ${UT_REQUIRED_LIBS})
add_test (ut_QuickGridHorizon python run_unit_test.py ./ut_QuickGridHorizon)

