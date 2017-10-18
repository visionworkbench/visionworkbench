file(GLOB VW_GEOMETRY_SRC *.h *.cc *.hpp *.tcc)
vw_add_library(vw_geometry "vw_core;vw_math;${Boost_LIBRARIES}" ${VW_GEOMETRY_SRC})
