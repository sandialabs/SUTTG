#
# CMakeLists.txt file for the SUTTG examples directory
#

set(exampleinput
waterfullopt97-1p0fast.xml
waterfullopt97-1p0opta.xml
waterfullopt97-1p0optb.xml
waterfullopt97-0p1optb.xml
waterfullopt95-1p0fast.xml
waterhighp95-0p01optb.xml
)

#
# install examples to share/suttg/examples
#
install(FILES ${exampleinput} README DESTINATION share/suttg/examples)
