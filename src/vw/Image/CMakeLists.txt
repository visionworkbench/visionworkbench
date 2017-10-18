file(GLOB VW_IMAGE_SRC *.h *.cc *.tcc)
vw_add_library(vw_image "vw_core;vw_math;${Boost_LIBRARIES}" ${VW_IMAGE_SRC})
target_link_libraries(vw_image PRIVATE ${OpenCV_LIBS})