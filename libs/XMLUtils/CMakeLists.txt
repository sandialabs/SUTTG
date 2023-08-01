#--------------------------------------------------------------------#
#  Copyright (2013) Sandia Corporation. Under the terms of Contract  #
#  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    #
#  retains certain rights in this software.                          #
#                                                                    #
#  For full terms of distribution see the accompanying LICENSE file. #
#  To distribute this file outside of SUTTG, substitue the full text #
#  of the LICENSE file for this notice.                              #
#--------------------------------------------------------------------#

#
# CMakeLists.txt file for the XMLUtils Expat wrapper
#
project(XMLUtils)

#
# Source files
#
set(xmlutil_sources
XMLAttributeList.cpp  XMLElement.cpp  XMLExpatParser.cpp
XMLParser.cpp  XMLUtils.cpp
)

#
# Library
#
add_library(XMLUtils ${xmlutil_sources})

#
# Need expat
#
include(FindEXPAT)

if (EXPAT_FOUND)
include_directories(${EXPAT_INCLUDE_DIRS})
target_link_libraries(XMLUtils ${EXPAT_LIBRARIES})
else()
message(SEND_ERROR "XMLUtils error: could not find EXPAT library")
endif()
