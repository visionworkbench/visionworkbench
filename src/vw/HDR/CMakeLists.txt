file(GLOB VW_HDR_SRC *.h *.cc *.tcc)
vw_add_library(vw_hdr "vw_core;vw_image;vw_math;vw_camera;${Boost_LIBRARIES}" ${VW_HDR_SRC})
