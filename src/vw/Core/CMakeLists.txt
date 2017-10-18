file(GLOB VW_CORE_SRC *.h *.cc *.tcc)
vw_add_library(vw_core "${Boost_LIBRARIES}" ${VW_CORE_SRC})
