file(GLOB VW_STEREO_SRC *.h *.hpp *.cc *.tcc)
vw_add_library(vw_stereo "vw_core;vw_image;vw_cartography;vw_camera;${Boost_LIBRARIES}" ${VW_STEREO_SRC})
