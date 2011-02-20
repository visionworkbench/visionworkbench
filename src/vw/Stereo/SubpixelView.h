// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_AFFINE_SUBPIXEL_VIEW__
#define __VW_STEREO_AFFINE_SUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Interpolation.h>

namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  template <class PreprocFilterT, class ImageT>
  class SubpixelView : public ImageViewBase<SubpixelView<PreprocFilterT, ImageT> > {

    ImageViewRef<PixelMask<Vector2f> > m_disparity_map;
    ImageT m_left_image, m_right_image;

    // General Settings
    Vector2i m_kernel_size;
    bool m_do_h_subpixel, m_do_v_subpixel;
    int32 m_which_subpixel;
    PreprocFilterT m_preproc_filter;
    bool m_verbose;

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<SubpixelView> pixel_accessor;

    template <class DisparityViewT>
    SubpixelView(DisparityViewT const& disparity_map,
                 ImageT const& left_image,
                 ImageT const& right_image,
                 int32 kern_width, int32 kern_height,
                 bool do_horizontal_subpixel,
                 bool do_vertical_subpixel,
                 int32 which_subpixel,
                 PreprocFilterT preproc_filter,
                 bool verbose) : m_disparity_map(disparity_map),
                                 m_left_image(left_image),
                                 m_right_image(right_image),
                                 m_kernel_size(Vector2i(kern_width,kern_height)),
                                 m_do_h_subpixel(do_horizontal_subpixel),
                                 m_do_v_subpixel(do_vertical_subpixel),
                                 m_which_subpixel(which_subpixel),
                                 m_preproc_filter(preproc_filter),
                                 m_verbose(verbose) {
      // Basic assertions
      VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                (left_image.impl().rows() == right_image.impl().rows()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()),
                ArgumentErr() << "SubpixelView::SubpixelView(): input image dimensions and/or disparity_map dimensions do not agree.\n");

      VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                ArgumentErr() << "SubpixelView::SubpixelView(): multi-channel, multi-plane images not supported.\n");
    }

    // Standard ImageView interface methods
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(float /*x*/, float /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "SubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {

      // Find the range of disparity values for this patch.
      BBox2i search_range = get_disparity_range(crop(m_disparity_map, bbox));

      // The area in the right image that we'll be searching is
      // determined by the bbox of the left image plus the search
      // range.
      BBox2i left_crop_bbox(bbox);
      BBox2i right_crop_bbox(bbox.min() + search_range.min(),
                             bbox.max() + search_range.max());

      // The correlator requires the images to be the same size. The
      // search bbox will always be larger than the given left image
      // bbox, so we just make the left bbox the same size as the
      // right bbox.
      left_crop_bbox.max() = left_crop_bbox.min() + right_crop_bbox.size();

      // Finally, we must adjust both bounding boxes to account for
      // the size of the kernel itself.
      right_crop_bbox.min() -= m_kernel_size;
      right_crop_bbox.max() += m_kernel_size;
      left_crop_bbox.min() -=  m_kernel_size;
      left_crop_bbox.max() +=  m_kernel_size;

      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<float> left_image_patch, right_image_patch;
      // parabola subpixel does the same preprocessing as the pyramid correlator
      left_image_patch = m_preproc_filter(crop(edge_extend(m_left_image,ZeroEdgeExtension()),
                                               left_crop_bbox));
      right_image_patch = m_preproc_filter(crop(edge_extend(m_right_image,ZeroEdgeExtension()),
                                                right_crop_bbox));
      ImageView<PixelMask<Vector2f> > disparity_map_patch =
        crop(edge_extend(m_disparity_map, ZeroEdgeExtension()),
             left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      PixelMask<Vector2f> disparity_patch_translation( search_range.min() );
      disparity_map_patch -= disparity_patch_translation;

      switch (m_which_subpixel){

      case 0 : // No Subpixel
        break;
      case 1 : // Parabola Subpixel
        subpixel_correlation_parabola(disparity_map_patch,
                                      left_image_patch,
                                      right_image_patch,
                                      m_kernel_size[0], m_kernel_size[1],
                                      m_do_h_subpixel, m_do_v_subpixel,
                                      m_verbose);
        break;
      case 2: // Bayes EM  Subpixel
        {
          const int32 pyramid_levels = 2;
          std::vector< ImageView<float> > l_patches, r_patches;
          std::vector<BBox2i> rois;
          ImageView<PixelMask<Vector2f> > d_subpatch;

          // I'd like for image subsampling to use a gaussian when
          // downsampling however it was introducing some edge effects
          // that I couldn't figure out within a reasonable time frame.
          for ( int32 i = 0; i < pyramid_levels; i++ ) {
            if ( i > 0 ) {
              // Building all other levels
              l_patches.push_back( subsample( l_patches.back(), 2 ) );
              r_patches.push_back( subsample( r_patches.back(), 2 ) );
              ImageView<PixelMask<Vector2f> > d_subpatch_buf =
                disparity_subsample( d_subpatch );
              d_subpatch = d_subpatch_buf;
              rois.push_back( rois.back()/2 );
            } else {
              // First step down from native resolution
              l_patches.push_back( subsample( left_image_patch, 2 ) );
              r_patches.push_back( subsample( right_image_patch, 2 ) );
              d_subpatch = disparity_subsample( disparity_map_patch );
              rois.push_back( BBox2i( m_kernel_size[0], m_kernel_size[1],
                                      bbox.width(), bbox.height() ) );
            }
          }

          for ( int32 i = pyramid_levels-1; i >= 0; i-- ) {
            subpixel_optimized_affine_2d_EM(d_subpatch,
                                            l_patches[i], r_patches[i],
                                            m_kernel_size[0], m_kernel_size[1],
                                            rois[i],
                                            m_do_h_subpixel, m_do_v_subpixel,
                                            m_verbose);
            BBox2i crop_bbox;
            if ( i > 0 )
              crop_bbox = BBox2i(0,0,l_patches[i-1].cols(),
                                 l_patches[i-1].rows());
            else
              crop_bbox = BBox2i(0,0,left_image_patch.cols(),
                                 left_image_patch.rows());
            ImageView<PixelMask<Vector2f> > d_subpatch_buf =
              crop(disparity_upsample(edge_extend(d_subpatch)), crop_bbox);
            d_subpatch = d_subpatch_buf;
          }

          disparity_map_patch = d_subpatch;

          // Perfrom final pass at native resolution
          subpixel_optimized_affine_2d_EM(disparity_map_patch,
                                          left_image_patch,
                                          right_image_patch,
                                          m_kernel_size[0], m_kernel_size[1],
                                          BBox2i(m_kernel_size[0],m_kernel_size[1],
                                                 bbox.width(), bbox.height()),
                                          m_do_h_subpixel, m_do_v_subpixel,
                                          m_verbose);
        }
        break;
      default:
        vw_throw(ArgumentErr() << "Unknown subpixel correlation type: " << m_which_subpixel << ".");
        break;
      }

      // Undo the above adjustment
      disparity_map_patch += disparity_patch_translation;

      // This may seem confusing, but we must crop here so that the
      // good pixel data is placed into the coordinates specified by
      // the bbox.  This allows rasterize to touch those pixels
      // using the coordinates inside the bbox.  The pixels outside
      // those coordinates are invalid, and they never get accessed.
      return crop(disparity_map_patch, BBox2i(m_kernel_size[0]-bbox.min()[0],
                                              m_kernel_size[1]-bbox.min()[1],
                                              m_left_image.cols(),
                                              m_left_image.rows() ));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };

  template <class PreprocFilterT, class ImageT, class DisparityT>
  SubpixelView<PreprocFilterT, ImageT>
  subpixel_refine( ImageViewBase<DisparityT> const& disparity_map,
                   ImageViewBase<ImageT> const& left_image,
                   ImageViewBase<ImageT> const& right_image,
                   int32 kern_width, int32 kern_height,
                   bool do_horizontal, bool do_vertical,
                   int32 which_subpixel, PreprocFilterT const& filter,
                   bool verbose = false ) {
    typedef SubpixelView<PreprocFilterT,ImageT> result_type;
    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), kern_width, kern_height,
                        do_horizontal, do_vertical, which_subpixel, filter,
                        verbose );
  }

  template <class PreprocFilterT, class ImageT, class DisparityT>
  SubpixelView<PreprocFilterT, ImageT>
  subpixel_refine( ImageViewBase<DisparityT> const& disparity_map,
                   ImageViewBase<ImageT> const& left_image,
                   ImageViewBase<ImageT> const& right_image,
                   Vector2i const& kernel_size,
                   bool do_horizontal, bool do_vertical,
                   int32 which_subpixel, PreprocFilterT const& filter,
                   bool verbose = false ) {
    typedef SubpixelView<PreprocFilterT,ImageT> result_type;
    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), kernel_size[0], kernel_size[1],
                        do_horizontal, do_vertical, which_subpixel, filter,
                        verbose );
  }

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
