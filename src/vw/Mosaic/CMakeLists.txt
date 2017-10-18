file(GLOB VW_MOSAIC_SRC *.h *.cc *.tcc *.hpp)
vw_add_library(vw_mosaic "vw_core;vw_image;vw_fileio;vw_cartography;${Boosst_LIBARRIES}" ${VW_MOSAIC_SRC})
