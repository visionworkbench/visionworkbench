// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_STEREO_PARABOLASUBPIXEL_VIEW__
#define __VW_STEREO_PARABOLASUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Stereo/PrefilterEnum.h>

// TODO(oalexan1): The integer disparity is in fact float. Convert all to float.
namespace vw { namespace stereo {

  class ParabolaSubpixelView : public ImageViewBase<ParabolaSubpixelView> {
    ImageViewRef<PixelMask<Vector2f>> m_disparity;
    ImageViewRef<PixelGray<float>> m_left_image;
    ImageViewRef<PixelGray<float>> m_right_image;
    Vector2i   m_kernel_size;

    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See PrefilterEnum.h for the types
    float m_prefilter_width;     ///< Preprocessing filter width

    Matrix<float,6,9> m_p_A_matrix;

    /// Compute the subpixel disparity for each input integer disparity
    /// - This function is written in a roundabout way to maximize the benefit from our
    ///   fast_box_sum() function.  Of course, we already performed all of these computations
    ///   back when we generated the integer disparity!!!
    ImageView<PixelMask<Vector2f>>
    evaluate(ImageView<PixelMask<Vector2i>> const& integer_disparity,    // Cropped input disp
             ImageViewRef<PixelGray<float>> const& left_filtered_image,  // uncropped left
             ImageViewRef<PixelGray<float>> const& right_filtered_image, // uncropped right
             BBox2i const& left_region,      ///< The ROI of the left and right images
             BBox2i const& right_region,     ///  to use in order to do our computations.
             BBox2i const& disparity_region, ///< The ROI in the entire output image
             BBox2i const& search_range ) const; /// Range of input disparity values in region

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<ParabolaSubpixelView> pixel_accessor;

    ParabolaSubpixelView(ImageViewRef<PixelMask<Vector2f>>    const& disparity,
                         ImageViewRef<PixelGray<float>>    const& left_image,
                         ImageViewRef<PixelGray<float>>    const& right_image,
                         PrefilterModeType prefilter_mode, float prefilter_width,
                         Vector2i const& kernel_size ) :
      m_disparity( disparity.impl() ), m_left_image( left_image.impl() ),
      m_right_image( right_image.impl() ), 
      m_kernel_size( kernel_size ),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width)
      {
      VW_ASSERT( m_disparity.cols() == m_left_image.cols() &&
                 m_disparity.rows() == m_left_image.rows(),
                 ArgumentErr() << "SubpixelView: Disparity image must match left image." );

      // We get a considerable speedup in our 2d subpixel correlation if
      // we go ahead and compute the pseudoinverse of the A matrix (where
      // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
      // for the range of x = [-1:1] and y = [-1:1].
      static float pinvA_data[] =/*
        { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
          1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
          1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
          -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
          -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
          -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 };*/
        {  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,   1.0/6, -1.0/3,  1.0/6,  // = a
         1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,   1.0/6,  1.0/6,  1.0/6,  // = b
         1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0,  -1.0/4,    0.0,  1.0/4,  // = c
        -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6,  -1.0/6,    0.0,  1.0/6,  // = d
        -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,    0.0,   1.0/6,  1.0/6,  1.0/6,  // = e
        -1.0/9,  2.0/9, -1.0/9,  2.0/9,   5.0/9, 2.0/9,  -1.0/9,  2.0/9, -1.0/9 };// = f
      m_p_A_matrix = Matrix<float,6,9>( pinvA_data );
    }

    inline int32 cols  () const { return m_disparity.cols(); }
    inline int32 rows  () const { return m_disparity.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator() ( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0 ) const {
      vw_throw( NoImplErr() << "SubpixelView:operator() has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section does the actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    prerasterize_type prerasterize( BBox2i const& bbox ) const;

    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize(prerasterize(bbox), dest, bbox );
    }
  }; // End class ParabolaSubpixelView

  ParabolaSubpixelView
  parabola_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity,
                    ImageViewRef<PixelGray<float>> const& left_image,
                    ImageViewRef<PixelGray<float>> const& right_image,
                    PrefilterModeType prefilter_mode, float prefilter_width,
                    Vector2i const& kernel_size);
  
}} // namespace vw::stereo

#endif // __VW_STEREO_PARABOLASUBPIXEL_VIEW__
