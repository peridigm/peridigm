#
# Main executable for Peridigm
#

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
    ${Trilinos_LIBRARIES}
    ${PDNEIGH_LIBS}
    ${Zoltan_LIBRARY}
    ${PARSER_LIBS}
    ${MESH_INPUT_LIBS}
    ${UTILITIES_LIBS}
    ${VTK_LIBS}
    ${MESH_OUTPUT_LIBS}
    ${REQUIRED_LIBS}
    ${PdMaterialUtilitiesLib}
    ${Peridigm_LIBRARY}
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
