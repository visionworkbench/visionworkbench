// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_CORRELATOR_VIEW__
#define __VW_STEREO_CORRELATOR_VIEW__

#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PyramidCorrelator.h>
#include <vw/Stereo/DisparityMap.h>

#include <ostream>

namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  template <class ImagePixelT, class MaskPixelT, class PreProcFuncT>
  class CorrelatorView : public ImageViewBase<CorrelatorView<ImagePixelT, MaskPixelT, PreProcFuncT> > {

    ImageViewRef<ImagePixelT> m_left_image, m_right_image;
    ImageViewRef<MaskPixelT> m_left_mask, m_right_mask;
    PreProcFuncT m_preproc_func;

    // Settings
    BBox2i m_search_range;
    Vector2i m_kernel_size;
    float m_cross_corr_threshold;
    float m_corr_score_threshold;
    int32 m_cost_blur;
    stereo::CorrelatorType m_correlator_type;
    std::string m_debug_prefix;
    bool m_do_pyramid_correlator;

    // Precalculated constants
    int32 m_num_pyramid_levels;
    Vector2i m_kernpad;         // Padding used around a render box

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<CorrelatorView> pixel_accessor;

    template <class ImageT, class MaskT>
    CorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
                   ImageViewBase<MaskT> const& left_mask, ImageViewBase<MaskT> const& right_mask,
                   PreProcFuncT const& preproc_func,
                   bool do_pyramid_correlator = true ) :
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
      m_preproc_func(preproc_func), m_do_pyramid_correlator(do_pyramid_correlator) {

        // Basic assertions
        VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                  (left_image.impl().rows() == right_image.impl().rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image dimensions do not agree.\n");

        VW_ASSERT((left_image.impl().cols() == left_mask.impl().cols()) &&
                  (left_image.impl().rows() == left_mask.impl().rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image and mask image dimensions do not agree.\n");

        VW_ASSERT((left_image.impl().cols() == right_mask.impl().cols()) &&
                  (left_image.impl().rows() == right_mask.impl().rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image and mask image dimensions do not agree.\n");

        VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                  (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): multi-channel, multi-plane images not supported.\n");

        // Set some sensible default values
        m_search_range = BBox2i(-50,-50,100,100);
        m_kernel_size = Vector2i(24,24);
        m_cross_corr_threshold = 2.0;
        m_corr_score_threshold = 1.3;
        m_cost_blur = 1;
        m_correlator_type = ABS_DIFF_CORRELATOR;

        // Calculating constants
        m_num_pyramid_levels = 4;
        if ( !m_do_pyramid_correlator )
          m_num_pyramid_levels = 1;
        m_kernpad = m_kernel_size*pow(2,m_num_pyramid_levels-1)/2;
      }

      // Basic accessor functions
      void set_search_range(BBox2i range) { m_search_range = range; }
      BBox2i search_range() const { return m_search_range; }

      void set_kernel_size(Vector2i size) {
        m_kernel_size = size;
        m_kernpad = m_kernel_size*pow(2,m_num_pyramid_levels-1)/2;
      }
      Vector2i kernel_size() const { return m_kernel_size; }

      void set_correlator_options(int32 cost_blur, stereo::CorrelatorType correlator_type) {
        m_cost_blur = cost_blur;
        m_correlator_type = correlator_type;
      }
      int32 cost_blur() const { return m_cost_blur; }
      stereo::CorrelatorType correlator_type() const { return m_correlator_type; }

      void set_cross_corr_threshold(float threshold) { m_cross_corr_threshold = threshold; }
      float cross_corr_threshold() const { return m_cross_corr_threshold; }

      void set_corr_score_threshold(float threshold) { m_corr_score_threshold = threshold; }
      float corr_score_threshold() const { return m_corr_score_threshold; }

      /// Turn on debugging output.  The debug_file_prefix string is
      /// used as a prefix for all debug image files.
      void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }

      // Standard ImageView interface methods
      inline int32 cols() const { return m_left_image.cols(); }
      inline int32 rows() const { return m_left_image.rows(); }
      inline int32 planes() const { return 1; }

      inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

      inline pixel_type operator()(double /*i*/, double /*j*/, int32 /*p*/ = 0) const {
        vw_throw(NoImplErr() << "CorrelatorView::operator()(double i, double j, int32 p) has not been implemented.");
        return pixel_type();
      }

      /// \cond INTERNAL
      typedef CropView<ImageView<pixel_type> > prerasterize_type;
      inline prerasterize_type prerasterize(BBox2i bbox) const {
        vw_out(DebugMessage, "stereo") << "CorrelatorView: rasterizing image block " << bbox << ".\n";

        // The area in the right image that we'll be searching is
        // determined by the bbox of the left image plus the search
        // range.
        BBox2i left_crop_bbox(bbox);
        BBox2i right_crop_bbox( bbox.min() + m_search_range.min(),
                                bbox.max() + m_search_range.max() );

        // The correlator requires the images to be the same size. The
        // search bbox will always be larger than the given left image
        // bbox, so we just make the left bbox the same size as the
        // right bbox.
        left_crop_bbox.max() = left_crop_bbox.min() +
          Vector2i(right_crop_bbox.width(), right_crop_bbox.height());

        // Finally, we must adjust both bounding boxes to account for
        // the size of the kernel itself.
        right_crop_bbox.min() -= m_kernpad;
        right_crop_bbox.max() += m_kernpad;
        left_crop_bbox.min() -= m_kernpad;
        left_crop_bbox.max() += m_kernpad;

        // Log some helpful debugging info
        vw_out(DebugMessage, "stereo") << "\t search_range:    "
                                       << m_search_range << std::endl;
        vw_out(DebugMessage, "stereo") << "\t left_crop_bbox:  "
                                       << left_crop_bbox << std::endl;
        vw_out(DebugMessage, "stereo") << "\t right_crop_bbox: "
                                       << right_crop_bbox << std::endl;

        ImageView<pixel_type> disparity_map;

        { // Force the scope for the crops

          // We crop the images to the expanded bounding box and edge
          // extend in case the new bbox extends past the image bounds.
          //
          // We're not making a copy as the pyramid correlator will make its own copy
          ImageViewRef<ImagePixelT> cropped_left_image =
            crop(edge_extend(m_left_image, ZeroEdgeExtension()), left_crop_bbox);
          ImageViewRef<ImagePixelT> cropped_right_image =
            crop(edge_extend(m_right_image, ZeroEdgeExtension()), right_crop_bbox);
          ImageViewRef<MaskPixelT> cropped_left_mask =
            crop(edge_extend(m_left_mask, ZeroEdgeExtension()), left_crop_bbox);
          ImageViewRef<MaskPixelT> cropped_right_mask =
            crop(edge_extend(m_right_mask, ZeroEdgeExtension()), right_crop_bbox);

          // We have all of the settings adjusted.  Now we just have to
          // run the correlator.
          if ( m_do_pyramid_correlator ) {
            PyramidCorrelator correlator(BBox2(0,0,m_search_range.width(),
                                               m_search_range.height()),
                                         Vector2i(m_kernel_size[0], m_kernel_size[1]),
                                         m_cross_corr_threshold, m_corr_score_threshold,
                                         m_cost_blur, m_correlator_type, m_num_pyramid_levels);

            // For debugging: this saves the disparity map at various
            // pyramid levels to disk.
            if (!m_debug_prefix.empty()) {
              std::ostringstream ostr;
              ostr << "-" << bbox.min().x() << "-" << bbox.max().x()
                   << "_" << bbox.min().y() << "-" << bbox.max().y() << "-";
              correlator.set_debug_mode(m_debug_prefix + ostr.str());
            }

            disparity_map = correlator( cropped_left_image, cropped_right_image,
                                        cropped_left_mask, cropped_right_mask,
                                        m_preproc_func);
          } else {
            OptimizedCorrelator correlator(BBox2(0,0,m_search_range.width(),
                                                 m_search_range.height()),
                                           m_kernel_size[0],
                                           m_cross_corr_threshold, m_corr_score_threshold,
                                           m_cost_blur, m_correlator_type );
            disparity_map = disparity_mask(correlator( cropped_left_image,
                                                       cropped_right_image,
                                                       m_preproc_func ),
                                           cropped_left_mask,
                                           cropped_right_mask );
          }
        } // Ending scope of crops

        // Adjust the disparities to be relative to the uncropped
        // image pixel locations
        for (int32 v = 0; v < disparity_map.rows(); ++v)
          for (int32 u = 0; u < disparity_map.cols(); ++u)
            if (is_valid(disparity_map(u,v)) )
              remove_mask(disparity_map(u,v)) += m_search_range.min();

        // This may seem confusing, but we must crop here so that the
        // good pixel data is placed into the coordinates specified by
        // the bbox.  This allows rasterize to touch those pixels
        // using the coordinates inside the bbox.  The pixels outside
        // those coordinates are invalid, and they never get accessed.
        return CropView<ImageView<pixel_type> > (disparity_map, BBox2i(m_kernpad[0]-bbox.min().x(),
                                                                       m_kernpad[1]-bbox.min().y(),
                                                                       bbox.width(), bbox.height()));
      }

      template <class DestT>
      inline void rasterize(DestT const& dest, BBox2i bbox) const {
        vw::rasterize(prerasterize(bbox), dest, bbox);
      }
      /// \endcond
  };

  /// Summarize a CorrelatorView object
  template <class ImagePixelT, class MaskPixelT, class PreProcFuncT>
  std::ostream& operator<<( std::ostream& os, CorrelatorView<ImagePixelT,MaskPixelT,PreProcFuncT> const& view ) {
    os << "------------------------- CorrelatorView ----------------------\n";
    os << "\tsearch range: " << view.search_range() << "\n";
    os << "\tkernel size : " << view.kernel_size() << "\n";
    os << "\txcorr thresh: " << view.cross_corr_threshold() << "\n";
    os << "\tcost blur: " << view.cost_blur() << "\n";
    os << "\tcorrelator type: " << view.correlator_type() << "\n";
    os << "\tcorrscore rejection thresh: " << view.corr_score_threshold() << "\n";
    os << "---------------------------------------------------------------\n";
    return os;
  }

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
