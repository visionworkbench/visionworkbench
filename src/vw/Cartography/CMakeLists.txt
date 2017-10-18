file(GLOB VW_CARTOGRAPHY_SRC *.h *.cc *.tcc)
vw_add_library(vw_cartography "vw_core;vw_camera;vw_math;${Boost_LIBRARIES}" ${VW_CARTOGRAPHY_SRC})
target_link_libraries(vw_cartography PRIVATE ${GDAL_LIBRARY} ${Proj_LIBRARIES})
