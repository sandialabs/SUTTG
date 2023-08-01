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
# CMakeLists.txt file for the Sandia Unstructured Triangle Table Generator
#
project(SUTTG C CXX Fortran)

#
# library include directories
#
include_directories(${eosinc})
include_directories(${xmlinc})
include_directories(${netcdfinc})
include_directories(${utriinc})

#
# library link directories
#
link_directories(${netcdflib})

#
# need MPI
#
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

#
# pthreads required, Mac needs to use named semaphores
#
find_package(Threads REQUIRED)
if (APPLE)
   add_definitions(-DSEM_NAMED)
endif()

#
# Source files
#
set(generatorSources
buildtables.C BMesh.C TriMesh.C RMesh.C std_ncstatus.c
InputParser.C TMesh.C Tables.C UTriTables.C MCMin.C
)

#
# generator executable
#
add_executable(utriTableBuilder ${generatorSources})
#
# generator uses MPI
#
if(MPI_COMPILE_FLAGS)
  set_target_properties(utriTableBuilder PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(utriTableBuilder PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#
# link dependencies of generator
#
target_link_libraries(utriTableBuilder sesml utri XMLUtils qd ${MPI_LIBRARIES})

#
# table viz executable
#
add_executable(utriTableDump tableviz.C std_ncstatus.c)
target_link_libraries(utriTableDump utri ${MPI_LIBRARIES})

#
# delaunay triangulation code tester
#
include_directories(${qdinc})
add_executable(dtbtest deltribound_test.C)
target_link_libraries(dtbtest qd)

#
# testing follows
#

#
# test support file
#
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/redirect.sh "exec $1 >& $2\n") 

#
# macro for delaunay tests
#
macro(makedtbtest number)
configure_file(dtbtest.in${number} dtbtest.in${number} COPYONLY)
configure_file(dtbtest.out${number} dtbtest.out${number} COPYONLY)
add_test(NAME dtbtest${number} COMMAND sh redirect.sh "./dtbtest dtbtest.in${number}" dtbtest${number}.out)
add_test(NAME dtbtest${number}compare COMMAND ${CMAKE_COMMAND} -E compare_files dtbtest${number}.out dtbtest.out${number})
set_property(TEST dtbtest${number}compare APPEND PROPERTY DEPENDS dtbtest${number})
endmacro(makedtbtest)

#
# delaunay tests
#
makedtbtest(1)
makedtbtest(2)
makedtbtest(3)
makedtbtest(4)
makedtbtest(5)
makedtbtest(6)
makedtbtest(7)
makedtbtest(8)
makedtbtest(9)
makedtbtest(10)
makedtbtest(11)
makedtbtest(12)
makedtbtest(13)

#
# install executables
#
install(TARGETS utriTableBuilder utriTableDump RUNTIME DESTINATION bin)

