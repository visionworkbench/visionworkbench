file(GLOB VW_INTERESTPOINT_SRC *.h *.cc *.tcc)
vw_add_library(vw_interestpoint "vw_core;vw_image;vw_fileio;${Boost_LIBRARIES}" ${VW_INTERESTPOINT_SRC})
target_link_libraries(vw_interestpoint PRIVATE ${OpenCV_LIBS})
