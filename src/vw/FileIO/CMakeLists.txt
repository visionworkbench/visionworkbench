file(GLOB VW_FILEIO_SRC *.h *.tcc *.cc)
list(REMOVE_ITEM VW_FILEIO_SRC 
    ${CMAKE_CURRENT_LIST_DIR}/DiskImageResourceOpenEXR.cc
    ${CMAKE_CURRENT_LIST_DIR}/DiskImageResourceHDF.cc
    ${CMAKE_CURRENT_LIST_DIR}/MemoryImageResourceOpenEXR.cc
)
vw_add_library(vw_fileio "vw_core;vw_image;vw_math;${Boost_LIBRARIES}" ${VW_FILEIO_SRC})
target_link_libraries(vw_fileio PRIVATE ${GDAL_LIBRARY} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES})