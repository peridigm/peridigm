set(utPeridigm_Compute_Force_SOURCES
   ./utPeridigm_Compute_Force.cpp
)
set(utPeridigm_Compute_Angular_Momentum_SOURCES
   ./utPeridigm_Compute_Angular_Momentum.cpp
)
set(utPeridigm_Compute_Linear_Momentum_SOURCES
   ./utPeridigm_Compute_Linear_Momentum.cpp
)
#set(utPeridigm_Compute_Energy_SOURCES
#   ./utPeridigm_Compute_Energy.cpp
#)

set(utPeridigm_Compute_Kinetic_Energy_SOURCES
   ./utPeridigm_Compute_Kinetic_Energy.cpp
)
set(utPeridigm_Compute_Strain_Energy_Density_SOURCES
   ./utPeridigm_Compute_Strain_Energy_Density.cpp
)
set(utPeridigm_Compute_Strain_Energy_SOURCES
   ./utPeridigm_Compute_Strain_Energy.cpp
)

add_executable(utPeridigm_Compute_Force ${utPeridigm_Compute_Force_SOURCES})
target_link_libraries(utPeridigm_Compute_Force
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)
add_executable(utPeridigm_Compute_Angular_Momentum ${utPeridigm_Compute_Angular_Momentum_SOURCES})
target_link_libraries(utPeridigm_Compute_Angular_Momentum
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)
add_executable(utPeridigm_Compute_Linear_Momentum ${utPeridigm_Compute_Linear_Momentum_SOURCES})
target_link_libraries(utPeridigm_Compute_Linear_Momentum
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)
#add_executable(utPeridigm_Compute_Energy ${utPeridigm_Compute_Energy_SOURCES})
#target_link_libraries(utPeridigm_Compute_Energy
#   ${Peridigm_LIBRARY}
#   ${Peridigm_LINK_LIBRARIES}
#   ${Boost_LIBRARIES}
#)
add_executable(utPeridigm_Compute_Kinetic_Energy ${utPeridigm_Compute_Kinetic_Energy_SOURCES})
target_link_libraries(utPeridigm_Compute_Kinetic_Energy
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)
add_executable(utPeridigm_Compute_Strain_Energy_Density ${utPeridigm_Compute_Strain_Energy_Density_SOURCES})
target_link_libraries(utPeridigm_Compute_Strain_Energy_Density
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)
add_executable(utPeridigm_Compute_Strain_Energy ${utPeridigm_Compute_Strain_Energy_SOURCES})
target_link_libraries(utPeridigm_Compute_Strain_Energy
   ${Peridigm_LIBRARY}
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)

add_test (utPeridigm_Compute_Force python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Force)
add_test (utPeridigm_Compute_Force_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Force)

#add_test (utPeridigm_Compute_Angular_Momentum python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Angular_Momentum)
#add_test (utPeridigm_Compute_Angular_Momentum_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Angular_Momentum)

add_test (utPeridigm_Compute_Linear_Momentum python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Linear_Momentum)
add_test (utPeridigm_Compute_Linear_Momentum_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Linear_Momentum)

#add_test (utPeridigm_Compute_Energy python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Energy)
#add_test (utPeridigm_Compute_Energy_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Energy)

add_test (utPeridigm_Compute_Kinetic_Energy python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Kinetic_Energy)
add_test (utPeridigm_Compute_Kinetic_Energy_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Kinetic_Energy)

add_test (utPeridigm_Compute_Strain_Energy_Density python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Strain_Energy_Density)
add_test (utPeridigm_Compute_Strain_Energy_Density_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Strain_Energy_Density)

add_test (utPeridigm_Compute_Strain_Energy python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Strain_Energy)
add_test (utPeridigm_Compute_Strain_Energy_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Strain_Energy)
