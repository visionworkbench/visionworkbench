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
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/PhaseSubpixelView.h>

//#include <boost/foreach.hpp>

namespace vw {
namespace stereo {

  enum PyramidSubpixelView_Algorithm {
    SUBPIXEL_LUCAS_KANADE = 0,
    SUBPIXEL_FAST_AFFINE  = 1,
    SUBPIXEL_BAYES_EM     = 2,
    SUBPIXEL_PHASE        = 3
  };

  /// An image view for performing image correlation using affine sub-pixel correlation.
  template <class ImageT1, class ImageT2, class ImageTD>
  class PyramidSubpixelView : public ImageViewBase<PyramidSubpixelView<ImageT1, ImageT2, ImageTD> > {

    ImageTD m_disparity_map;
    ImageT1 m_left_image;
    ImageT2 m_right_image;

    // General Settings - Could any of these be shared with the other classes?
    Vector2i       m_kernel_size;
    int32          m_max_pyramid_levels;
    PyramidSubpixelView_Algorithm m_algorithm;

    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See Prefilter.h for the types
    float m_prefilter_width;     ///< Preprocessing filter width


  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type          result_type;
    typedef ProceduralPixelAccessor<PyramidSubpixelView> pixel_accessor;

    PyramidSubpixelView(ImageViewBase<ImageTD> const& disparity_map,
                        ImageViewBase<ImageT1> const& left_image,
                        ImageViewBase<ImageT2> const& right_image,
                        PrefilterModeType prefilter_mode, float prefilter_width,
                        Vector2i const& kernel_size,
                        int32 max_pyramid_levels,
                        PyramidSubpixelView_Algorithm algorithm) :
      m_disparity_map(disparity_map.impl()),
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_kernel_size(kernel_size), 
      m_max_pyramid_levels(max_pyramid_levels), m_algorithm(algorithm),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width) {

      // Basic assertions
      VW_ASSERT( m_disparity_map.cols() == m_left_image.cols() &&
                 m_disparity_map.rows() == m_left_image.rows(),
                 ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView():  Disparity image must match left image.\n");

      VW_ASSERT((m_left_image.channels() == 1) && (m_left_image.planes() == 1) &&
                (m_right_image.channels() == 1) && (m_right_image.planes() == 1),
                ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView(): multi-channel, multi-plane images are not supported.\n");

      // Max pyramid can't go below 0 ... or bayes em won't process anything
      if ( m_max_pyramid_levels < 0 ) m_max_pyramid_levels = 0;
    }

    // Standard ImageView interface methods
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(float /*x*/, float /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "PyramidSubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif
      // Find the range of disparity values for this patch.
      ImageView<pixel_type > disparity_map_patch =
        crop(m_disparity_map, bbox);
      BBox2i search_range = get_disparity_range(disparity_map_patch);

      // The area in the right image that we'll be searching is
      // determined by the bbox of the left image plus the search range.
      BBox2i left_crop_bbox(bbox);
      BBox2i right_crop_bbox(bbox.min() + search_range.min(),
                             bbox.max() + search_range.max());

      // TODO: Where in the code is this required?
      // The correlator requires the images to be the same size. The
      // search bbox will always be larger than the given left image
      // bbox, so we just make the left bbox the same size as the right bbox.
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
      // Using the convenient prefilter function yields slightly different results here
      //  for some reason so we use this method to avoid having to redo all the tests.
      if (m_prefilter_mode == PREFILTER_LOG) {
        stereo::LaplacianOfGaussian prefilter(m_prefilter_width);
        left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
        right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
      } else {
        if (m_prefilter_mode == PREFILTER_MEANSUB) {
          stereo::SubtractedMean prefilter(m_prefilter_width);
          left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
          right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
        } else { // PREFILTER_NONE
          stereo::NullOperation prefilter;
          left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
          right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
        }
      }
      disparity_map_patch = crop(edge_extend(m_disparity_map, ZeroEdgeExtension()),
                                 left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      PixelMask<Vector2f> disparity_patch_translation( search_range.min() );
      disparity_map_patch -= disparity_patch_translation;

      std::vector< ImageView<float> > l_patches, r_patches;
      std::vector<BBox2i> rois;
      ImageView<pixel_type > d_subpatch;

      // This initialization makes a difference only if the number of
      // levels is 0, which is a valid situation.
      d_subpatch = disparity_map_patch;

      // In general we should blur our images when downsampling, but when tested in
      //  this code it produced a small benefit in the unit tests but no benefit in any
      //  of the real images it was tested with.  Because of this it was deemed not
      //  worth the computation time to include it here.
      
      BBox2i full_res_roi( m_kernel_size[0], m_kernel_size[1],
                           bbox.width(), bbox.height() );
      
      // Create an image pyramid for the left and right image patches we will correlate.
      for ( int32 i = 0; i < m_max_pyramid_levels; i++ ) {
        if ( i > 0 ) {
          // Building all other levels
          l_patches.push_back( subsample( l_patches.back(), 2 ) );
          r_patches.push_back( subsample( r_patches.back(), 2 ) );
          ImageView<PixelMask<Vector2f> > d_subpatch_buf = disparity_subsample( d_subpatch );
          d_subpatch = d_subpatch_buf;
          rois.push_back( rois.back()/2 );
        } else {
          // First step down from native resolution
          l_patches.push_back( subsample( left_image_patch,  2 ) );
          r_patches.push_back( subsample( right_image_patch, 2 ) );
          d_subpatch = disparity_subsample( disparity_map_patch );
          rois.push_back(full_res_roi/2);
        }
      }

      const int PHASE_SUBPIXEL_ACCURACY = 20; // Theoretically accurate to 1/N of a pixel.

      // Loop through all but final pyramid levels.
      for ( int32 i = m_max_pyramid_levels-1; i >= 0; i-- ) {

        switch(m_algorithm) {
        case SUBPIXEL_LUCAS_KANADE:
          subpixel_optimized_LK_2d(d_subpatch,
                                   l_patches[i], r_patches[i],
                                   m_kernel_size[0], m_kernel_size[1],
                                   rois[i], true, true, false );
          break;
        case SUBPIXEL_FAST_AFFINE:
          subpixel_optimized_affine_2d(d_subpatch,
                                       l_patches[i], r_patches[i],
                                       m_kernel_size[0], m_kernel_size[1],
                                       rois[i], true, true, false );
          break;
        case SUBPIXEL_BAYES_EM:
          subpixel_optimized_affine_2d_EM(d_subpatch,
                                          l_patches[i], r_patches[i],
                                          m_kernel_size[0], m_kernel_size[1],
                                          rois[i], true, true, false );
          break;
        case SUBPIXEL_PHASE:
          subpixel_phase_2d(d_subpatch,
                            l_patches[i], r_patches[i],
                            m_kernel_size[0], m_kernel_size[1],
                            rois[i], PHASE_SUBPIXEL_ACCURACY);
          //write_image("/home/smcmich1/data/subpixel/disp_subpix_x.tif",
          //    select_channel(d_subpatch,0));
          break;
        default:
          vw_throw(NoImplErr() << "Invalid algorithm selection passed to PyramidSubpixelView.");
        }

        BBox2i crop_bbox;
        if ( i > 0 )
          crop_bbox = BBox2i(0,0,l_patches[i-1].cols(), l_patches[i-1].rows());
        else
          crop_bbox = BBox2i(0,0,left_image_patch.cols(), left_image_patch.rows());
        ImageView<pixel_type > d_subpatch_buf =
          crop(disparity_upsample(edge_extend(d_subpatch)), crop_bbox);
        d_subpatch = d_subpatch_buf;
      } // End loop through pyramid levels

      disparity_map_patch = d_subpatch;

      // Perfrom final pass at native resolution
      switch(m_algorithm) {
      case SUBPIXEL_LUCAS_KANADE:
        subpixel_optimized_LK_2d(disparity_map_patch,
                                 left_image_patch, right_image_patch,
                                 m_kernel_size[0], m_kernel_size[1],
                                 full_res_roi,
                                 true, true, false );
        break;
      case SUBPIXEL_FAST_AFFINE:
        subpixel_optimized_affine_2d(disparity_map_patch,
                                     left_image_patch, right_image_patch,
                                     m_kernel_size[0], m_kernel_size[1],
                                     full_res_roi,
                                     true, true, false );
        break;
      case SUBPIXEL_BAYES_EM:
        subpixel_optimized_affine_2d_EM(disparity_map_patch,
                                        left_image_patch, right_image_patch,
                                        m_kernel_size[0], m_kernel_size[1],
                                        full_res_roi,
                                        true, true, false );
        break;
      case SUBPIXEL_PHASE:
        subpixel_phase_2d(disparity_map_patch,
                         left_image_patch, right_image_patch,
                         m_kernel_size[0], m_kernel_size[1],
                         full_res_roi,
                         PHASE_SUBPIXEL_ACCURACY);
        break;
      default:
        vw_throw(NoImplErr() << "Invalid algorithm selection passed to PyramidSubpixelView.");
      }

      // Undo the above adjustment
      disparity_map_patch += disparity_patch_translation;

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << " s\n";
#endif

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

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };

  // Set of wrapper functions to help use PyramidSubpixelView
  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  lk_subpixel( ImageViewBase<DisparityT> const& disparity_map,
               ImageViewBase<ImageT1   > const& left_image,
               ImageViewBase<ImageT2   > const& right_image,
               PrefilterModeType prefilter_mode, float prefilter_width,
               Vector2i const& kernel_size,
               int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), 
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_LUCAS_KANADE);
  }

  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  affine_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                   ImageViewBase<ImageT1   > const& left_image,
                   ImageViewBase<ImageT2   > const& right_image,
                   PrefilterModeType prefilter_mode, float prefilter_width,
                   Vector2i const& kernel_size,
                   int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), 
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_FAST_AFFINE);
  }

  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  bayes_em_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                     ImageViewBase<ImageT1   > const& left_image,
                     ImageViewBase<ImageT2   > const& right_image,
                     PrefilterModeType prefilter_mode, float prefilter_width,
                     Vector2i const& kernel_size,
                     int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(),
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_BAYES_EM);
  }

  // Phase subpixel seems to work better without multi-resolution, so
  //  the default number of pyramid levels is zero.
  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  phase_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                  ImageViewBase<ImageT1   > const& left_image,
                  ImageViewBase<ImageT2   > const& right_image,
                  PrefilterModeType prefilter_mode, float prefilter_width,
                  Vector2i const& kernel_size,
                  int max_pyramid_levels = 0 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(),
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_PHASE);
  }

  // End of components for Pyramid subpixel view









}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
