macro(vw_add_library VW_LIB_NAME VW_DEPENDS)
    add_library(${VW_LIB_NAME} ${ARGN})
    target_link_libraries(${VW_LIB_NAME} PUBLIC ${VW_DEPENDS})
    install(
        TARGETS ${VW_LIB_NAME}
        EXPORT "vw-targets"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        RUNTIME DESTINATION "bin"
        INCLUDES DESTINATION "include"
    )
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DESTINATION include/vw
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.tcc"
        )
    set_target_properties(${VW_LIB_NAME} PROPERTIES FOLDER "vw-modules")
endmacro(vw_add_library)

if(VW_HAVE_PKG_CORE)
    add_subdirectory(Core)
endif(VW_HAVE_PKG_CORE)

if(VW_HAVE_PKG_MATH)
    add_subdirectory(Math)
endif(VW_HAVE_PKG_MATH)

if(VW_HAVE_PKG_IMAGE)
    add_subdirectory(Image)
endif(VW_HAVE_PKG_IMAGE)

if(VW_HAVE_PKG_FILEIO)
    add_subdirectory(FileIO)
endif(VW_HAVE_PKG_FILEIO)

if(VW_HAVE_PKG_INTERESTPOINT)
    add_subdirectory(InterestPoint)
endif(VW_HAVE_PKG_INTERESTPOINT)

if(VW_HAVE_PKG_CAMERA)
    add_subdirectory(Camera)
endif(VW_HAVE_PKG_CAMERA)

if(VW_HAVE_PKG_CARTOGRAPHY)
    add_subdirectory(Cartography)
endif(VW_HAVE_PKG_CARTOGRAPHY)

if(VW_HAVE_PKG_STEREO)
    add_subdirectory(Stereo)
endif(VW_HAVE_PKG_STEREO)

if(VW_HAVE_PKG_BUNDLEADJUSTMENT)
    add_subdirectory(BundleAdjustment)
endif(VW_HAVE_PKG_BUNDLEADJUSTMENT)

if(VW_HAVE_PKG_HDR)
    add_subdirectory(HDR)
endif(VW_HAVE_PKG_HDR)

if(VW_HAVE_PKG_MOSAIC)
    add_subdirectory(Mosaic)
endif(VW_HAVE_PKG_MOSAIC)

if(VW_HAVE_PKG_GEOMETRY)
    add_subdirectory(Geometry)
endif(VW_HAVE_PKG_GEOMETRY)

if(VW_HAVE_PKG_TOOL)
    add_subdirectory(tools)
endif(VW_HAVE_PKG_TOOL)
