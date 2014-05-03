
# --- Bayer application --------------------------------------------



#message(${FULL_LIBRARY_LIST})

# ---------------------------------------------------------------
# Define each of the tools to build

#TODO: Turn them off when components are not present

# Applies an inverse Bayer filter to a raw camera image
add_executable(bayer bayer.cc) 
target_link_libraries(bayer ${FULL_LIBRARY_LIST}) 

# Blends multiple images into a composite
add_executable(blend blend.cc) 
target_link_libraries(blend ${FULL_LIBRARY_LIST}) 

# Converts a DEM to a simple color map
add_executable(colormap colormap.cc Common.h Common.cc) 
target_link_libraries(colormap ${FULL_LIBRARY_LIST}) 

# Contour generation program
add_executable(contourgen contour.cpp FitCurves.cpp contourgen.cc contour.h) 
target_link_libraries(contourgen ${FULL_LIBRARY_LIST}) 

# Apply the block stereo correlator to two images, producing a disparity map
add_executable(correlate correlate.cc) 
target_link_libraries(correlate ${FULL_LIBRARY_LIST}) 

# TODO: This is built with multiple PIXEL_TYPE and CHANNEL_TYPE inputs!
# Creates a blended composite from multiple georeferenced images
#add_executable(geoblend geoblend.cc geoblend.h geoblend_help.cc Common.h Common.cc) 
#target_link_libraries(geoblend ${FULL_LIBRARY_LIST}) 

# Adds or adjusts an image's georeferencing information
add_executable(georef georef.cc) 
target_link_libraries(georef ${FULL_LIBRARY_LIST}) 

# A GPU-based version of the Stereo correlator
#add_executable(gpu_correlate gpu_correlate.cc) 
#target_link_libraries(gpu_correlate ${FULL_LIBRARY_LIST}) 

# Fades an image's alpha channel for blending
add_executable(grassfirealpha grassfirealpha.cc) 
target_link_libraries(grassfirealpha ${FULL_LIBRARY_LIST}) 

# Merge multiple source images into an HDR image
add_executable(hdr_merge hdr_merge.cc) 
target_link_libraries(hdr_merge ${FULL_LIBRARY_LIST}) 

# Map an HDR image into an 8-bit dynamic range
add_executable(hdr_tonemap hdr_tonemap.cc) 
target_link_libraries(hdr_tonemap ${FULL_LIBRARY_LIST}) 

# Converts a DEM to a simple shaded-relief map
add_executable(hillshade hillshade.cc Common.h Common.cc) 
target_link_libraries(hillshade ${FULL_LIBRARY_LIST}) 

#TODO: Replicate multi-library thing?
# Generates a quadtree of some sort from a georeferenced image
#add_executable(image2qtree image2qtree.cc image2qtree_help.cc) 
#target_link_libraries(image2qtree ${FULL_LIBRARY_LIST}) 

# Creates a TOAST quadtree from multiple georeferenced images
add_executable(image2toast image2toast.cc) 
target_link_libraries(image2toast ${FULL_LIBRARY_LIST}) 

add_executable(ipalign ipalign.cc) 
target_link_libraries(ipalign ${FULL_LIBRARY_LIST}) 

add_executable(ipfind ipfind.cc) 
target_link_libraries(ipfind ${FULL_LIBRARY_LIST}) 

add_executable(ipmatch ipmatch.cc) 
target_link_libraries(ipmatch ${FULL_LIBRARY_LIST}) 

add_executable(learnpca learnpca.cc) 
target_link_libraries(learnpca ${FULL_LIBRARY_LIST}) 

add_executable(make_ba_test_data make_ba_test_data.cc) 
target_link_libraries(make_ba_test_data ${FULL_LIBRARY_LIST}) 

# Prints out an image's EXIF headers
add_executable(print_exif print_exif.cc) 
target_link_libraries(print_exif ${FULL_LIBRARY_LIST}) 

# Calculates aspect and gradient angle for DEM
add_executable(slopemap slopemap.cc Common.h Common.cc) 
target_link_libraries(slopemap ${FULL_LIBRARY_LIST}) 








