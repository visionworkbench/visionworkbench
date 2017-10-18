file(GLOB VW_MATH_SRC *.h *.cc *.tcc)
vw_add_library(vw_math "vw_core;${Boost_LIBRARIES}" ${VW_MATH_SRC})
target_link_libraries(vw_math PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${Flann_LIBRARIES})
