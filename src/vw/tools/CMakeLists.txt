macro(vw_add_executable VW_EXE_NAME VW_EXE_MODULE)
    add_executable(${VW_EXE_NAME} ${ARGN})
    target_link_libraries(${VW_EXE_NAME} ${VW_EXE_MODULE})
    set_target_properties(${VW_EXE_NAME} PROPERTIES FOLDER "vw_tools")
    install(
        TARGETS ${VW_EXE_NAME}
        RUNTIME DESTINATION "bin"
    )
endmacro(vw_add_executable)

if(VW_HAVE_PKG_CAMERA)
    # undistort_image
    vw_add_executable(undistort_image vw_cartography undistort_image.cc)
    vw_add_executable(convert_pinhole_model vw_camera convert_pinhole_model.cc)
endif(VW_HAVE_PKG_CAMERA)

if(VW_HAVE_PKG_CARTOGRAPHY)
    vw_add_executable(hillshade vw_cartography hillshade.cc hillshade.h)
    vw_add_executable(georef vw_cartography georef.cc)
    vw_add_executable(colormap vw_cartography colormap.cc)
    vw_add_executable(georefquery vw_cartography georef_query.cc)
    vw_add_executable(slopemap vw_cartography slopemap.cc)
    vw_add_executable(grassfirealpha vw_cartography grassfirealpha.cc)
endif(VW_HAVE_PKG_CARTOGRAPHY)

if(VW_HAVE_PKG_HDR)
    vw_add_executable(hdr_tonemap vw_hdr hdr_tonemap.cc)
    vw_add_executable(hdr_merge vw_hdr hdr_merge.cc)
endif(VW_HAVE_PKG_HDR)

if(VW_HAVE_PKG_INTERESTPOINT AND VW_HAVE_PKG_LAPACK)
    vw_add_executable(ipalign vw_interestpoint ipalign.cc)
    vw_add_executable(ipmatch vw_interestpoint ipmatch.cc)
    vw_add_executable(ipfind vw_interestpoint ipfind.cc)
    vw_add_executable(learnpca vw_interestpoint learnpca.cc)
endif(VW_HAVE_PKG_INTERESTPOINT AND VW_HAVE_PKG_LAPACK)

# TODO get image2qtree done, this requires generate multiple object from a single template cc file

if(VW_HAVE_PKG_STEREO AND VW_HAVE_PKG_INTERESTPOINT)
    vw_add_executable(correlate "vw_interestpoint;vw_stereo" correlate.cc)
endif(VW_HAVE_PKG_STEREO AND VW_HAVE_PKG_INTERESTPOINT)

if(VW_HAVE_PKG_FLOOD_DETECT)
    vw_add_executable(detect_water vw_cartography detect_water.cc multispectral.h radar.h flood_common.h)
    vw_add_executable(clean_sentinel1_borders vw_cartography clean_sentinel1_borders.cc flood_common.h)
endif(VW_HAVE_PKG_FLOOD_DETECT)
