// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
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
  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreProcFuncT>
  class CorrelatorView : public ImageViewBase<CorrelatorView<Image1T,Image2T,Mask1T,Mask2T,PreProcFuncT> > {

    Image1T m_left_image;
    Image2T m_right_image;
    Mask1T m_left_mask;
    Mask2T m_right_mask;
    PreProcFuncT m_preproc_func;

    // Settings
    BBox2i m_search_range;
    Vector2i m_kernel_size;
    float m_cross_corr_threshold;
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

    CorrelatorView(ImageViewBase<Image1T> const& left_image,
                   ImageViewBase<Image2T> const& right_image,
                   ImageViewBase<Mask1T> const& left_mask,
                   ImageViewBase<Mask2T> const& right_mask,
                   PreProcFuncT const& preproc_func,
                   bool do_pyramid_correlator = true ) :
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
      m_preproc_func(preproc_func), m_do_pyramid_correlator(do_pyramid_correlator) {

        // Basic assertions
        VW_ASSERT((m_left_image.cols() == m_right_image.cols()) &&
                  (m_left_image.rows() == m_right_image.rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image dimensions do not agree.\n");

        VW_ASSERT((m_left_image.cols() == m_left_mask.cols()) &&
                  (m_left_image.rows() == m_left_mask.rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image and mask image dimensions do not agree.\n");

        VW_ASSERT((m_left_image.cols() == m_right_mask.cols()) &&
                  (m_left_image.rows() == m_right_mask.rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image and mask image dimensions do not agree.\n");

        VW_ASSERT((m_left_image.channels() == 1) && (m_left_image.planes() == 1) &&
                  (m_right_image.channels() == 1) && (m_right_image.planes() == 1),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): multi-channel, multi-plane images not supported.\n");

        // Set some sensible default values
        m_search_range = BBox2i(-50,-50,100,100);
        m_kernel_size = Vector2i(23,23);
        m_cross_corr_threshold = 2.0;
        m_cost_blur = 1;
        m_correlator_type = ABS_DIFF_CORRELATOR;

        // Calculating constants
        m_num_pyramid_levels = 4;
        if ( !m_do_pyramid_correlator )
          m_num_pyramid_levels = 1;
        m_kernpad = m_kernel_size*pow(2,m_num_pyramid_levels-1)/2;
      }

      // Basic accessor functions
      void set_search_range(BBox2i const& range) { m_search_range = range; }
      BBox2i search_range() const { return m_search_range; }

      void set_kernel_size(Vector2i const& size) {
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
      inline prerasterize_type prerasterize(BBox2i const& bbox) const {
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
        VW_OUT(DebugMessage, "stereo") << "\t search_range:    "
                                       << m_search_range << std::endl
                                       << "\t left_crop_bbox:  "
                                       << left_crop_bbox << std::endl
                                       << "\t right_crop_bbox: "
                                       << right_crop_bbox << std::endl;

        // We crop the images to the expanded bounding box and edge
        // extend in case the new bbox extends past the image bounds.
        ImageView<typename Mask1T::pixel_type> cropped_left_mask =
          crop(edge_extend(m_left_mask, ZeroEdgeExtension()), left_crop_bbox);
        ImageView<typename Mask2T::pixel_type> cropped_right_mask =
          crop(edge_extend(m_right_mask, ZeroEdgeExtension()), right_crop_bbox);

        // The result that we return
        ImageView<pixel_type> disparity_map(cropped_left_mask.cols(),
                                            cropped_left_mask.rows());

        if ( sum_of_pixel_values(cropped_left_mask) != 0 &&
             sum_of_pixel_values(cropped_right_mask) != 0 ) {
          // We have all of the settings adjusted.  Now we just have to
          // run the correlator.
          if ( m_do_pyramid_correlator ) {
            PyramidCorrelator correlator(BBox2i(0,0,m_search_range.width(),
                                                m_search_range.height()),
                                         m_kernel_size,
                                         m_cross_corr_threshold,
                                         m_cost_blur, m_correlator_type, m_num_pyramid_levels);

            // For debugging: this saves the disparity map at various
            // pyramid levels to disk.
            if (!m_debug_prefix.empty()) {
              std::ostringstream ostr;
              ostr << "-" << bbox.min().x() << "-" << bbox.max().x()
                   << "_" << bbox.min().y() << "-" << bbox.max().y() << "-";
              correlator.set_debug_mode(m_debug_prefix + ostr.str());
            }

            disparity_map =
              correlator( crop(edge_extend(m_left_image, ZeroEdgeExtension()), left_crop_bbox),
                          crop(edge_extend(m_right_image, ZeroEdgeExtension()), right_crop_bbox),
                          cropped_left_mask, cropped_right_mask,
                          m_preproc_func);
          } else {
            OptimizedCorrelator correlator(BBox2i(0,0,m_search_range.width(),
                                                  m_search_range.height()),
                                           m_kernel_size[0],
                                           m_cross_corr_threshold,
                                           m_cost_blur, m_correlator_type );
            disparity_map =
              disparity_mask(correlator( crop(edge_extend(m_left_image, ZeroEdgeExtension()), left_crop_bbox),
                                         crop(edge_extend(m_right_image, ZeroEdgeExtension()), right_crop_bbox),
                                         m_preproc_func ),
                             cropped_left_mask,
                             cropped_right_mask );
          }
        } // Ending scope of crops

        // Adjust the disparities to be relative to the uncropped
        // image pixel locations
        // This should just be a straight forward add
        disparity_map += pixel_type(m_search_range.min());

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
      inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
        vw::rasterize(prerasterize(bbox), dest, bbox);
      }
      /// \endcond
  };

  /// Helper function so one doesn't have to type the expanded type
  template <class View1T, class View2T, class Mask1T, class Mask2T, class PreProcT>
  CorrelatorView<View1T,View2T,Mask1T,Mask2T,PreProcT>
  correlate( ImageViewBase<View1T> const& left_image,
             ImageViewBase<View2T> const& right_image,
             ImageViewBase<Mask1T> const& left_mask,
             ImageViewBase<Mask2T> const& right_mask,
             PreProcT const& preproc_func,
             BBox2i const& search_range = BBox2i(-50,-50,100,100),
             Vector2i const& kernel = Vector2i(23,23),
             stereo::CorrelatorType correlator_type = stereo::ABS_DIFF_CORRELATOR,
             int32 cost_blur = 1,
             bool do_pyramid_correlator = true ) {
    typedef CorrelatorView<View1T,View2T,Mask1T,Mask2T,PreProcT> result_type;
    result_type result( left_image.impl(), right_image.impl(),
                        left_mask.impl(), right_mask.impl(),
                        preproc_func, do_pyramid_correlator );
    result.set_search_range( search_range );
    result.set_kernel_size( kernel );
    result.set_correlator_options(cost_blur, correlator_type);
    return result;
  }


  /// Summarize a CorrelatorView object
  template <class Image1T, class Image2T, class Mask1T, class Mask2T, class PreProcFuncT>
  std::ostream& operator<<( std::ostream& os, CorrelatorView<Image1T,Image2T,Mask1T,Mask2T,PreProcFuncT> const& view ) {
    os << "------------------------- CorrelatorView ----------------------\n";
    os << "\tsearch range: " << view.search_range() << "\n";
    os << "\tkernel size : " << view.kernel_size() << "\n";
    os << "\txcorr thresh: " << view.cross_corr_threshold() << "\n";
    os << "\tcost blur: " << view.cost_blur() << "\n";
    os << "\tcorrelator type: " << view.correlator_type() << "\n";
    os << "---------------------------------------------------------------\n";
    return os;
  }

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
