file(GLOB VW_BUNDLEADJUSTMENT_SRC *.h *.cc *.tcc *.hpp)
vw_add_library(vw_bundleadjustment "vw_core;vw_math;vw_cartography;vw_fileio;vw_stereo;vw_interestpoint" ${VW_BUNDLEADJUSTMENT_SRC})
