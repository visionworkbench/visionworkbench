

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

# Adds or adjusts an image's georeferencing information
add_executable(georef georef.cc) 
target_link_libraries(georef ${COMMON_LIBS}) 
install(TARGETS georef DESTINATION bin)

# Fades an image's alpha channel for blending
add_executable(grassfirealpha grassfirealpha.cc) 
target_link_libraries(grassfirealpha ${COMMON_LIBS}) 
install(TARGETS grassfirealpha DESTINATION bin)

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

# image2qtree - Generates a quadtree of some sort from a georeferenced image
# - It takes forever to compile due to templating
#   so we split it up into multiple targets for parallelization.
add_library(libimage2qtree_grayau8  STATIC image2qtree_help.cc)
add_library(libimage2qtree_grayau16 STATIC image2qtree_help.cc)
add_library(libimage2qtree_grayai16 STATIC image2qtree_help.cc)
add_library(libimage2qtree_grayaf32 STATIC image2qtree_help.cc)
add_library(libimage2qtree_rgbau8   STATIC image2qtree_help.cc)
add_library(libimage2qtree_rgbau16  STATIC image2qtree_help.cc)
add_library(libimage2qtree_rgbai16  STATIC image2qtree_help.cc)
add_library(libimage2qtree_rgbaf32  STATIC image2qtree_help.cc)
target_compile_definitions(libimage2qtree_grayau8  PRIVATE -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=uint8  )
target_compile_definitions(libimage2qtree_grayau16 PRIVATE -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=uint16 )
target_compile_definitions(libimage2qtree_grayai16 PRIVATE -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=int16  )
target_compile_definitions(libimage2qtree_grayaf32 PRIVATE -DPIXEL_TYPE=PixelGrayA -DCHANNEL_TYPE=float32)
target_compile_definitions(libimage2qtree_rgbau8   PRIVATE -DPIXEL_TYPE=PixelRGBA  -DCHANNEL_TYPE=uint8  )
target_compile_definitions(libimage2qtree_rgbau16  PRIVATE -DPIXEL_TYPE=PixelRGBA  -DCHANNEL_TYPE=uint16 )
target_compile_definitions(libimage2qtree_rgbai16  PRIVATE -DPIXEL_TYPE=PixelRGBA  -DCHANNEL_TYPE=int16  )
target_compile_definitions(libimage2qtree_rgbaf32  PRIVATE -DPIXEL_TYPE=PixelRGBA  -DCHANNEL_TYPE=float32)
# Now that we have defined the individual static library files, link the executable to all of them.
add_executable(image2qtree image2qtree.cc)
target_link_libraries(image2qtree ${COMMON_LIBS} VwMosaic
                      libimage2qtree_grayau8  libimage2qtree_grayau16
                      libimage2qtree_grayai16 libimage2qtree_grayaf32
                      libimage2qtree_rgbau8   libimage2qtree_rgbau16
                      libimage2qtree_rgbai16  libimage2qtree_rgbaf32 )
install(TARGETS image2qtree DESTINATION bin)



