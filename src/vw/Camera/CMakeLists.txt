file(GLOB VW_CAMERA_SRC *.h *.cc *.tcc)
vw_add_library(vw_camera "vw_core;vw_math;vw_interestpoint;${Boost_LIBRARIES}" ${VW_CAMERA_SRC})
