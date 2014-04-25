
# --- Bayer application --------------------------------------------

# TODO: Don't just use this for all tools
set(FULL_LIBRARY_LIST
	${LAPACK_LIBRARIES} 
	${JPEG_LIBRARIES} 
	${GDAL_LIBRARY}
	${PROJ4_LIBRARY}
	${OPENEXR_LIBRARIES}
    VwCore
    VwImage
    VwFileIO
    VwCamera
	${Boost_LIBRARIES} 
)

#message(${FULL_LIBRARY_LIST})

add_executable(bayer bayer.cc) 

target_link_libraries(bayer ${FULL_LIBRARY_LIST}) 









