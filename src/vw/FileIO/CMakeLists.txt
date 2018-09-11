

# List files manually for this complicated folder
 
 
if (VW_HAVE_PKG_PNG)
  SET(png_headers DiskImageResourcePNG.h  MemoryImageResourcePNG.h  PngIO.h)
  SET(png_sources DiskImageResourcePNG.cc MemoryImageResourcePNG.cc PngIO.cc)
endif()

if (VW_HAVE_PKG_JPEG)
  SET(jpeg_headers DiskImageResourceJPEG.h  MemoryImageResourceJPEG.h  JpegIO.h)
  SET(jpeg_sources DiskImageResourceJPEG.cc MemoryImageResourceJPEG.cc JpegIO.cc)
endif()

if (VW_HAVE_PKG_TIFF)
  SET(tiff_headers DiskImageResourceTIFF.h)
  SET(tiff_sources DiskImageResourceTIFF.cc)
endif()

if (VW_HAVE_PKG_OPENEXR)
  if (ENABLE_EXCEPTIONS)
    SET(openexr_headers DiskImageResourceOpenEXR.h MemoryImageResourceOpenEXR.h)
    SET(openexr_sources DiskImageResourceOpenEXR.cc MemoryImageResourceOpenEXR.cc)
  endif()
endif()

if (VW_HAVE_PKG_HDF)
  SET(hdf_headers DiskImageResourceHDF.h)
  SET(hdf_sources DiskImageResourceHDF.cc)
endif()

if (VW_HAVE_PKG_GDAL)
  SET(gdal_headers DiskImageResourceGDAL.h  MemoryImageResourceGDAL.h  GdalIO.h)
  SET(gdal_sources DiskImageResourceGDAL.cc MemoryImageResourceGDAL.cc GdalIO.cc)
endif()
 
 
set(VW_FILEIO_SRC_FILES 
    DiskImageManager.h
    DiskImageResource.h
    DiskImageResource_internal.h
    DiskImageResourcePBM.h 
    DiskImageResourcePDS.h 
    DiskImageResourceRaw.h
    DiskImageUtils.h 
    DiskImageView.h 
    FileUtils.h
    FileUtils.cc
    MemoryImageResource.h 
    KML.h 
    ScanlineIO.h 
    TemporaryFile.h 
    ${gdal_headers} 
    ${hdf_headers} 
    ${jpeg_headers} 
    ${openexr_headers} 
    ${png_headers} 
    ${tiff_headers}
    DiskImageResource.cc 
    DiskImageResourcePBM.cc 
    DiskImageResourcePDS.cc 
    DiskImageResourceRaw.cc
    KML.cc 
    MemoryImageResource.cc 
    ScanlineIO.cc 
    TemporaryFile.cc 
    ${gdal_sources} 
    ${hdf_sources} 
    ${jpeg_sources} 
    ${openexr_sources} 
    ${png_sources} 
    ${tiff_sources}
    )


# Use wrapper function at this level to avoid code duplication
add_library_wrapper(VwFileIO "${VW_FILEIO_SRC_FILES}" "${VW_FILEIO_TEST_FILES}" "${VW_FILEIO_LIB_DEPENDENCIES}")





