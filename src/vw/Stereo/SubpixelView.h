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


#ifndef __VW_STEREO_SUBPIXEL_VIEW__
#define __VW_STEREO_SUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/PrefilterEnum.h>

namespace vw { namespace stereo {

  enum PyramidSubpixelView_Algorithm {
    SUBPIXEL_LUCAS_KANADE = 0,
    SUBPIXEL_FAST_AFFINE  = 1,
    SUBPIXEL_BAYES_EM     = 2,
    SUBPIXEL_PHASE        = 3
  };

  /// An image view for performing image correlation using affine sub-pixel correlation.
  class PyramidSubpixelView : public ImageViewBase<PyramidSubpixelView> {

    ImageViewRef<PixelMask<Vector2f>> m_disparity_map;
    ImageViewRef<PixelGray<float>>    m_left_image;
    ImageViewRef<PixelGray<float>>    m_right_image;

    // General Settings - Could any of these be shared with the other classes?
    Vector2i       m_kernel_size;
    int32          m_max_pyramid_levels;
    PyramidSubpixelView_Algorithm m_algorithm;

    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See PrefilterEnum.h for the types
    float m_prefilter_width;            ///< Preprocessing filter width
    int m_phase_subpixel_accuracy;


  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type          result_type;
    typedef ProceduralPixelAccessor<PyramidSubpixelView> pixel_accessor;
    PyramidSubpixelView(ImageViewRef<PixelMask<Vector2f>> const& disparity_map,
                        ImageViewRef<PixelGray<float>>    const& left_image,
                        ImageViewRef<PixelGray<float>>    const& right_image,
                        PrefilterModeType prefilter_mode, float prefilter_width,
                        Vector2i const& kernel_size,
                        int32 max_pyramid_levels,
                        PyramidSubpixelView_Algorithm algorithm,
                        int32 phase_subpixel_accuracy = 20) :
      m_disparity_map(disparity_map.impl()),
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_kernel_size(kernel_size), 
      m_max_pyramid_levels(max_pyramid_levels), m_algorithm(algorithm),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width),
      m_phase_subpixel_accuracy(phase_subpixel_accuracy) {

      // Basic assertions
      VW_ASSERT(m_disparity_map.cols() == m_left_image.cols() &&
                m_disparity_map.rows() == m_left_image.rows(),
                ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView(): "
                << "Disparity image must match left image.\n");

      VW_ASSERT((m_left_image.channels() == 1) && (m_left_image.planes() == 1) &&
                (m_right_image.channels() == 1) && (m_right_image.planes() == 1),
                ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView(): "
                << "multi-channel, multi-plane images are not supported.\n");

      // Max pyramid can't go below 0 ... or bayes em won't process anything
      if (m_max_pyramid_levels < 0) m_max_pyramid_levels = 0;
    }

    // Standard ImageView interface methods
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }

    inline pixel_type operator()(float /*x*/, float /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "PyramidSubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    prerasterize_type prerasterize(BBox2i const& bbox) const;

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };
  
  // Set of wrapper functions to help use PyramidSubpixelView
  PyramidSubpixelView
  lk_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity_map,
              ImageViewRef<PixelGray<float>> const& left_image,
              ImageViewRef<PixelGray<float>> const& right_image,
              PrefilterModeType prefilter_mode, float prefilter_width,
              Vector2i const& kernel_size,
              int max_pyramid_levels = 2);
  
  PyramidSubpixelView
  affine_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity_map,
                  ImageViewRef<PixelGray<float>> const& left_image,
                  ImageViewRef<PixelGray<float>> const& right_image,
                  PrefilterModeType prefilter_mode, float prefilter_width,
                  Vector2i const& kernel_size,
                  int max_pyramid_levels = 2);

  PyramidSubpixelView
  bayes_em_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity_map,
                    ImageViewRef<PixelGray<float>> const& left_image,
                    ImageViewRef<PixelGray<float>> const& right_image,
                    PrefilterModeType prefilter_mode, float prefilter_width,
                    Vector2i const& kernel_size,
                    int max_pyramid_levels = 2);
  
  // Phase subpixel seems to work better without multi-resolution, so
  //  the default number of pyramid levels is zero.
  PyramidSubpixelView
  phase_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity_map,
                 ImageViewRef<PixelGray<float>> const& left_image,
                 ImageViewRef<PixelGray<float>> const& right_image,
                 PrefilterModeType prefilter_mode, float prefilter_width,
                 Vector2i const& kernel_size,
                 int max_pyramid_levels = 0,
                 int phase_subpixel_accuracy = 20);

  // End of components for Pyramid subpixel view

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
