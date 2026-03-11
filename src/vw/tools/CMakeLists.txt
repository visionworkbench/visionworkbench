

# ---------------------------------------------------------------
# Define each of the tools to build


# A few header files from this folder are used by ASP.
INSTALL(FILES "Common.h"    DESTINATION include/vw/tools)

set(COMMON_LIBS VwCore VwMath VwFileIO VwImage VwCartography)

#TODO: Turn them off when components are not present

#TODO: Add flood detection components

# Converts a DEM to a simple color map
add_executable(colormap colormap.cc Common.h Common.cc) 
target_link_libraries(colormap ${COMMON_LIBS}) 
install(TARGETS colormap DESTINATION bin)

# Convert a pinhole model from one type to another
add_executable(convert_pinhole_model convert_pinhole_model.cc) 
target_link_libraries(convert_pinhole_model ${COMMON_LIBS} VwCamera) 
install(TARGETS convert_pinhole_model DESTINATION bin)

# Apply the block stereo correlator to two images, producing a disparity map
add_executable(correlate correlate.cc) 
target_link_libraries(correlate ${COMMON_LIBS} VwStereo VwInterestPoint) 
install(TARGETS correlate DESTINATION bin)

# Converts a DEM to a simple shaded-relief map
add_executable(hillshade hillshade.cc Common.h Common.cc)
target_link_libraries(hillshade ${COMMON_LIBS})
install(TARGETS hillshade DESTINATION bin)

add_executable(ipfind ipfind.cc)
target_link_libraries(ipfind ${COMMON_LIBS} VwInterestPoint)
install(TARGETS ipfind DESTINATION bin)

add_executable(ipmatch ipmatch.cc)
target_link_libraries(ipmatch ${COMMON_LIBS} VwInterestPoint)
install(TARGETS ipmatch DESTINATION bin)

# Prints out an image's EXIF headers
add_executable(print_exif print_exif.cc)
target_link_libraries(print_exif ${COMMON_LIBS})
install(TARGETS print_exif DESTINATION bin)

# Calculates aspect and gradient angle for DEM
add_executable(slopemap slopemap.cc Common.h Common.cc)
target_link_libraries(slopemap ${COMMON_LIBS})
install(TARGETS slopemap DESTINATION bin)

# Undistorts a pinhole image.
add_executable(undistort_image undistort_image.cc)
target_link_libraries(undistort_image ${COMMON_LIBS})
install(TARGETS undistort_image DESTINATION bin)

# image2qtree - Generates a quadtree from a georeferenced image
add_executable(image2qtree image2qtree.cc)
target_link_libraries(image2qtree ${COMMON_LIBS} VwMosaic)
install(TARGETS image2qtree DESTINATION bin)



