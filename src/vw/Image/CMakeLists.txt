 
# This is a complicated folder so we need to select files manually

if (HAVE_PKG_OPENCV)
  SET(opencv_headers ImageResourceOpenCV.h)
  SET(opencv_sources ImageResourceOpenCV.cc)
endif()
 
 
set(VW_IMAGE_SRC_FILES 
                      AlgorithmFunctions.h 
                      Algorithms.h 
                      BlockProcessor.h 
                      BlockRasterize.h 
                      Convolution.h 
                      EdgeExtend.h 
                      EdgeExtension.h 
                      Filter.h 
                      Filter.tcc 
                      ImageIO.h 
                      ImageMath.h 
                      ImageResource.h 
                      ImageResourceImpl.h 
                      ImageResourceStream.h 
                      ImageResourceView.h 
                      ImageViewBase.h 
                      ImageView.h 
                      ImageViewRef.h 
                      Interpolation.h 
                      Manipulation.h 
                      MaskViews.h 
                      Palette.h 
                      PerPixelAccessorViews.h 
                      PerPixelViews.h 
                      PixelAccessors.h 
                      PixelIterator.h 
                      PixelMask.h 
                      PixelMath.h 
                      PixelTypeInfo.h 
                      PixelTypes.h 
                      SparseImageCheck.h 
                      Statistics.h 
                      Transform.h 
                      UtilityViews.h 
                      ViewImageResource.h
                      ${opencv_headers} 
                      Filter.cc 
                      ImageResource.cc 
                      ImageResourceStream.cc 
                      Interpolation.cc 
                      PixelTypeInfo.cc
                      ${opencv_sources} 
     )


# Use wrapper function at this level to avoid code duplication
add_library_wrapper(VwImage "${VW_IMAGE_SRC_FILES}" "${VW_IMAGE_TEST_FILES}" "${VW_IMAGE_LIB_DEPENDENCIES}")


