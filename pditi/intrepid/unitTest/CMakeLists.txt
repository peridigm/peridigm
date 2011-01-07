#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#


# Add test
add_executable(utLinkIntrepid utLinkIntrepid.cxx)
target_link_libraries(utLinkIntrepid ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utLinkIntrepid python run_unit_test.py ./utLinkIntrepid)

add_executable(utFIFE_2D utFIFE_2D.cxx)
target_link_libraries(utFIFE_2D IIntrepid ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utFIFE_2D python run_unit_test.py ./utFIFE_2D)

add_executable(utTeuchosSmartPointer utTeuchosSmartPointer.cxx)
target_link_libraries(utTeuchosSmartPointer IIntrepid ${PERIDIGM_LIBS} ${VTK_LIBS}  ${UT_REQUIRED_LIBS} ${TRILINOS_LIBS})
add_test (utTeuchosSmartPointer python run_unit_test.py ./utTeuchosSmartPointer)

