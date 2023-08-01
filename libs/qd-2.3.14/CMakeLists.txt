#
# CMakeLists.txt file for the qd library
#

#
# Source files
#
set( qd_sources

src/dd_real.cpp src/dd_const.cpp src/fpu.cpp src/qd_real.cpp
src/qd_const.cpp src/util.cpp src/bits.cpp

)

#
# Library
#
add_library(qd ${qd_sources})
include_directories(${qdinc})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
