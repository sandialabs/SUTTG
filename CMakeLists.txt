cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release CACHE STRING "")

#
# CMakeLists.txt file for the SUTTG project
#
project(SUTTG C CXX Fortran)

enable_testing()

#
# define include directories
# before the add_subdirectory calls so they are inherited in scope
#
set(eosinc ${PROJECT_SOURCE_DIR}/libs/eos)
set(qdinc ${PROJECT_SOURCE_DIR}/libs/qd-2.3.14/include)
set(xmlinc ${PROJECT_SOURCE_DIR}/libs/XMLUtils)
set(utriinc ${PROJECT_SOURCE_DIR}/libs/utri ${PROJECT_BINARY_DIR}/libs/utri)
set(gtestinc ${PROJECT_SOURCE_DIR}/libs/gtest-1.6.0/include)

#
# save the base directory
#
set(suttg_source_dir ${PROJECT_SOURCE_DIR})

#
# site dependent include locations
#
include(config/config.cmake)

#
# compiler flags
#
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "-Wall -Wextra" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O2" CACHE STRING "" FORCE)
endif()
if(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_C_FLAGS "-Wall -Wextra" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_RELEASE "-O2" CACHE STRING "" FORCE)
endif()
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-Wall -Wextra" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING "" FORCE)
endif()

#
# add desired project directories
#
add_subdirectory(libs)
add_subdirectory(uteos)
add_subdirectory(examples)
add_subdirectory(data)

#
# install text files
#
install(FILES LICENSE README DESTINATION share/suttg)

#
# Source packaging
#
set(CPACK_SOURCE_GENERATOR "ZIP;TGZ")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "suttg-0.1")

set(CPACK_SOURCE_IGNORE_FILES ".*~" "[.]svn" "buildopt" "tests")

include(CPack)
