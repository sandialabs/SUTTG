#
# CMakeLists.txt file for the SUTTG data directory
#

set(datafiles
h2o95-1p0fast-pe.nc
h2o97-0p1optb-pe.nc
h2o97-1p0fast-pe.nc
h2o97-1p0opta-pe.nc
h2o97-1p0optb-pe.nc
h2o95highp-0p01optb-pe.nc
)

#
# install data to share/suttg/data
#
install(FILES ${datafiles} DESTINATION share/suttg/data)
