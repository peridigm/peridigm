#
# Main executable for Peridigm
#

# Glob header files in subdirectories
file (GLOB computeHEADERS 
  "${COMPUTE_DIR}/*.hpp"
)
file (GLOB contactHEADERS 
  "${CONTACT_DIR}/*.hpp"
)
file (GLOB coreHEADERS 
  "${CORE_DIR}/*.hpp"
)
file (GLOB discretizationHEADERS 
  "${DISCRETIZATION_DIR}/*.hpp"
)
file (GLOB damageHEADERS 
  "${DAMAGE_DIR}/*.hpp"
)
file (GLOB ioHEADERS 
  "${IO_DIR}/*.hpp"
)
file (GLOB materialsHEADERS 
  "${MATERIALS_DIR}/*.hpp"
)
file (GLOB quickGridHEADERS 
  "${QUICK_GRID_DIR}/*.h"
)
file (GLOB utilitiesHEADERS 
  "${UTILITIES_DIR}/*.h*"
)
file (GLOB pdneighHEADERS 
  "${PDNEIGH_DIR}/*.h"
)
file (GLOB parserHEADERS
  "${PARSER_DIR}/*.h"
)

set(Peridigm_HEADERS
    ${COMPUTE_DIR}/compute_includes.hpp
    ${computeHEADERS}
    ${contactHEADERS}
    ${coreHEADERS}
    ${discretizationHEADERS}
    ${damageHEADERS}
    ${ioHEADERS}
    ${materialsHEADERS}
    ${quickGridHEADERS}
    ${utilitiesHEADERS}
    ${pdneighHEADERS}
    ${parserHEADERS}
)

# Glob source files in subdirectories
file (GLOB computeSOURCES 
  "${COMPUTE_DIR}/*.cpp"
)
file (GLOB contactSOURCES 
  "${CONTACT_DIR}/*.cpp"
)
file (GLOB coreSOURCES 
  "${CORE_DIR}/*.cpp"
)
file (GLOB discretizationSOURCES 
  "${DISCRETIZATION_DIR}/*.cpp"
)
file (GLOB damageSOURCES 
  "${DAMAGE_DIR}/*.cpp"
)
file (GLOB ioSOURCES 
  "${IO_DIR}/*.cpp"
)
file (GLOB materialsSOURCES 
  "${MATERIALS_DIR}/*.cpp"
)

set(Peridigm_SOURCES
    ${COMPUTE_DIR}/compute_includes.hpp
    ${computeSOURCES}
    ${contactSOURCES}
    ${coreSOURCES}
    ${discretizationSOURCES}
    ${damageSOURCES}
    ${ioSOURCES}
    ${materialsSOURCES}
)

add_library(PeridigmLib ${Peridigm_SOURCES})
set(Peridigm_LIBRARY PeridigmLib)

set(Peridigm_LINK_LIBRARIES
    ${LCM_LIBRARY}
    ${Peridigm_LIBRARY}
    ${Trilinos_LIBRARIES}
    ${PDNEIGH_LIBS}
    ${Zoltan_LIBRARY}
    ${PARSER_LIBS}
    ${MESH_INPUT_LIBS}
    ${UTILITIES_LIBS}
    ${REQUIRED_LIBS}
    ${PdMaterialUtilitiesLib}
    ${Boost_LIBRARIES}
)

add_executable(Peridigm ${CORE_DIR}/Peridigm_Main.cpp)
target_link_libraries(Peridigm ${Peridigm_LINK_LIBRARIES})

#
# Install target for Peridigm main executable
#
install(TARGETS Peridigm
    RUNTIME DESTINATION bin
)
set_property(TARGET Peridigm
    PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE
)

# 
# Add subdirectories
#
add_subdirectory (compute/)
add_subdirectory (core/)
add_subdirectory (evaluators/)
add_subdirectory (io/)
add_subdirectory (materials/)
add_subdirectory (muParser/)

IF (INSTALL_PERIDIGM)
  configure_package_config_file(PeridigmConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PeridigmConfig.cmake
    INSTALL_DESTINATION ${LIB_INSTALL_DIR}/Peridigm/cmake
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
  install(FILES ${Peridigm_HEADERS}
    DESTINATION ${INCLUDE_INSTALL_DIR})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PeridigmConfig.cmake
    DESTINATION ${LIB_INSTALL_DIR}/Peridigm/cmake )
  install(TARGETS PeridigmLib EXPORT peridigm-export
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}) 
  FILE(WRITE "${CMAKE_BINARY_DIR}/peridigm-targets.cmake"
    "set(Peridigm_LIBRARIES ${Peridigm_LIBRARY} ${PDNEIGH_LIBS} ${PARSER_LIBS} ${MESH_INPUT_LIBS} ${UTILITIES_LIBS} ${PdMaterialUtilitiesLib})")
  FILE(INSTALL "${CMAKE_BINARY_DIR}/peridigm-targets.cmake"
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}/Peridigm/cmake)
ENDIF ()
