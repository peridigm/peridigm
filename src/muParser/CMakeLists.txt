#
# File for muParser library
#

# muParser library
add_library(muParser muParserBase.cpp  muParserBytecode.cpp muParserCallback.cpp  muParser.cpp  muParserError.cpp  muParserTokenReader.cpp)

IF (INSTALL_PERIDIGM)
   install(TARGETS muParser EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()