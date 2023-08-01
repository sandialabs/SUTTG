#
# CMakeLists.txt file for the SUTTG project libs directory
#

#
# help sesml find gtest
#
set(SESML_GTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gtest-1.6.0/include)
set(SESML_GTEST_IN_TREE 1)

#
# add desired libraries
#
add_subdirectory(eos)
add_subdirectory(gtest-1.6.0)
add_subdirectory(qd-2.3.14)
add_subdirectory(utri)
add_subdirectory(XMLUtils)
