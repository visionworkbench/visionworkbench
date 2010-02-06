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
  template <class PreprocFilterT>
  class SubpixelView : public ImageViewBase<SubpixelView<PreprocFilterT> > {

    ImageViewRef<PixelMask<Vector2f> > m_disparity_map;
    ImageViewRef<float> m_left_image;
    ImageViewRef<float> m_right_image;

    // General Settings
    int m_kern_width, m_kern_height;
    bool m_do_h_subpixel, m_do_v_subpixel;
    //bool m_do_affine_subpixel;
    int m_do_affine_subpixel;
    PreprocFilterT m_preproc_filter;
    bool m_verbose;

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<SubpixelView> pixel_accessor;
      

    //image subsampling by two
    ImageView<float> subsample_img_by_two(ImageView<float> &img) const {

      //determine the new size of the image
      int new_width = img.cols()/2;
      int new_height = img.rows()/2;

      ImageView<float> outImg(new_width, new_height, img.planes());
      ImageViewRef<float> interpImg = interpolate(img);

      for (vw::int32 p = 0; p < outImg.planes(); p++)
        for (vw::int32 i = 0; i < outImg.cols(); i++)
          for (vw::int32 j = 0; j < outImg.rows(); j++)
            outImg(i,j,p) = interpImg(2.0*i, 2.0*j, p);

      #if 0
      ImageView<float> g_img;
      g_img = gaussian_filter(img, 1.5);
      vw_out() << "Gaussian blurring\n";

      for (vw::int32 p = 0; p < outImg.planes(); p++)
        for (vw::int32 i = 0; i < outImg.cols(); i++)
          for (vw::int32 j = 0; j < outImg.rows(); j++)
            outImg(i,j,p) = g_img(2*i, 2*j, p);

      #endif

      return outImg;
    }

    // disparity map down-sampling by two
    ImageView<PixelMask<Vector2f> >
    subsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp) const {

      // determine the new size of the image
      int new_width = input_disp.cols()/2;
      int new_height = input_disp.rows()/2;

      ImageView<PixelMask<Vector2f> > outDisp(new_width, new_height);
      ImageViewRef<PixelMask<Vector2f> > disp = edge_extend(input_disp, ConstantEdgeExtension() );

      for (vw::int32 j = 0; j < outDisp.rows(); j++) {
        for (vw::int32 i = 0; i < outDisp.cols(); i++) {

          vw::int32 old_i = i*2;
          vw::int32 old_j = j*2;

          PixelMask<Vector2f> sum(0,0);
          int num_valid = 0;

          if ( is_valid(disp(old_i,old_j)) ) {
            sum += disp(old_i,old_j);
            num_valid++;
          }
          if ( is_valid(disp(old_i,old_j+1)) ) {
            sum += disp(old_i,old_j+1);
            num_valid++;
          }
          if ( is_valid(disp(old_i+1,old_j)) ) {
            sum += disp(old_i+1,old_j);
            num_valid++;
          }
          if ( is_valid(disp(old_i+1,old_j+1)) ) {
            sum += disp(old_i+1,old_j+1);
            num_valid++;
          }

          if (num_valid == 0)
            invalidate( outDisp(i,j) );
          else
            outDisp(i,j) = sum/float(2*num_valid);
        }
      }

      return outDisp;
    }

    // disparity map up-sampling by two
    ImageView<PixelMask<Vector2f> >
    upsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp,
                             int up_width, int up_height) const {

      ImageView<PixelMask<Vector2f> > outDisp(up_width, up_height);
      ImageViewRef<PixelMask<Vector2f> > disp = edge_extend(input_disp, ConstantEdgeExtension());

      for (vw::int32 j = 0; j < outDisp.rows(); ++j) {
        for (vw::int32 i = 0; i < outDisp.cols(); ++i) {
          float x = math::impl::_floor(float(i)/2.0), y = math::impl::_floor(float(j)/2.0);

          if ( i%2 == 0 && j%2 == 0) {
            if ( !is_valid(disp(x,y)) )
              invalidate( outDisp(i,j) );
            else
              outDisp(i,j) = 2 * disp(x,y);
          }

          else if (j%2 == 0) {
            if ( !is_valid(disp(x,y)) && !is_valid(disp(x+1,y)) )
              invalidate(outDisp(i,j));
            else {
              if ( is_valid(disp(x,y)) && is_valid(disp(x+1,y)) ) {
                outDisp(i,j) = disp(x,y) + disp(x+1,y);
              } else if ( is_valid(disp(x,y)) && !is_valid(disp(x+1,y)) ) {
                outDisp(i,j) = 2 * disp(x,y);
              } else if ( !is_valid(disp(x,y)) && is_valid(disp(x+1,y)) ) {
                outDisp(i,j) = 2 * disp(x+1,y);
              }
            }
          }

          else if (i%2 == 0) {
            if ( !is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) )
              invalidate( outDisp(i,j) );
            else {
              if ( is_valid(disp(x,y)) && is_valid(disp(x,y+1)) ) {
                outDisp(i,j) = disp(x,y) + disp(x,y+1);
              } else if ( is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) ) {
                outDisp(i,j) = 2*disp(x,y);
              } else if ( !is_valid(disp(x,y)) && is_valid(disp(x,y+1)) ) {
                outDisp(i,j) = 2*disp(x,y+1);
              }
            }
          }

          else {
            if ( is_valid(disp(x,y)) && is_valid(disp(x,y+1)) &&
                 is_valid(disp(x+1,y)) && is_valid(disp(x+1,y+1)) ) {

              // All good pixels
              float normx = float(i)/2.0-x, normy = float(j)/2.0-y, norm1mx = 1.0-normx, norm1my = 1.0-normy;
              outDisp(i,j) = PixelMask<Vector2f>( 2 * (disp(x  ,y  )[0] * norm1mx*norm1my +
                                                       disp(x+1,y  )[0] * normx*norm1my +
                                                       disp(x  ,y+1)[0] * norm1mx*normy +
                                                       disp(x+1,y+1)[0] * normx*normy),
                                                  2 * (disp(x  ,y  )[1] * norm1mx*norm1my +
                                                       disp(x+1,y  )[1] * normx*norm1my +
                                                       disp(x  ,y+1)[1] * norm1mx*normy +
                                                       disp(x+1,y+1)[1] * normx*normy) );

            } else if ( !is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) &&
                        !is_valid(disp(x+1,y)) && !is_valid(disp(x+1,y+1)) ) {
              // no good pixels
              invalidate( outDisp(i,j) );

            } else {
              // some good & some bad pixels
              //
              // We fudge things a little bit here by picking the
              // first good pixel.  This isn't exactly correct, but
              // it's close enough in the handful of cases where we
              // are near some missing pixels, and we just need an
              // approximately valid value.  The subpixel refinement
              // will correct any minor mistakes we introduce here.
              if ( is_valid(disp(x,y)) ) {
                outDisp(i,j) = 2*disp(x,y);
              } else if ( is_valid(disp(x+1,y)) ) {
                outDisp(i,j) = 2*disp(x+1,y);
              } else if ( is_valid(disp(x,y+1)) ) {
                outDisp(i,j) = 2*disp(x,y+1);
              } else if ( is_valid(disp(x+1,y+1)) ) {
                outDisp(i,j) = 2*disp(x+1,y+1);
              }
            }
          }
        }
      }

      return outDisp;
    }

    template <class DisparityViewT, class InputViewT>
    SubpixelView(DisparityViewT const& disparity_map,
                 InputViewT const& left_image,
                 InputViewT const& right_image,
                 int kern_width, int kern_height,
                 bool do_horizontal_subpixel,
                 bool do_vertical_subpixel,
                 //bool do_affine_subpixel,
                 int do_affine_subpixel,
                 PreprocFilterT preproc_filter,
                 bool verbose) : m_disparity_map(disparity_map),
                                 m_left_image(left_image),
                                 m_right_image(right_image),
                                 m_kern_width(kern_width), m_kern_height(kern_height),
                                 m_do_h_subpixel(do_horizontal_subpixel),
                                 m_do_v_subpixel(do_vertical_subpixel),
                                 m_do_affine_subpixel(do_affine_subpixel),
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

    inline pixel_type operator()(float x, float y, int32 p = 0) const {
      vw_throw(NoImplErr() << "SubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {

      // Find the range of disparity values for this patch.
      BBox2i search_range;
      try {
        search_range = get_disparity_range(crop(m_disparity_map, bbox));
      } catch ( std::exception &e ) {
        search_range = BBox2i();
      }

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
      left_crop_bbox.max() = left_crop_bbox.min() + Vector2i(right_crop_bbox.width(), right_crop_bbox.height());

      // Finally, we must adjust both bounding boxes to account for
      // the size of the kernel itself.
      right_crop_bbox.min() -= Vector2i(m_kern_width, m_kern_height);
      right_crop_bbox.max() += Vector2i(m_kern_width, m_kern_height);
      left_crop_bbox.min() -= Vector2i(m_kern_width, m_kern_height);
      left_crop_bbox.max() += Vector2i(m_kern_width, m_kern_height);

      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<float> left_image_patch, right_image_patch;
      if (m_do_affine_subpixel > 1) {
        // affine subpixel does its own pre-processing
        left_image_patch = crop(edge_extend(m_left_image,ZeroEdgeExtension()),
                                left_crop_bbox);
        right_image_patch = crop(edge_extend(m_right_image,ZeroEdgeExtension()),
                                 right_crop_bbox);
      } else {
        // parabola subpixel does the same preprocessing as the pyramid correlator
        left_image_patch = m_preproc_filter(crop(edge_extend(m_left_image,ZeroEdgeExtension()),
                                                 left_crop_bbox));
        right_image_patch = m_preproc_filter(crop(edge_extend(m_right_image,ZeroEdgeExtension()),
                                                  right_crop_bbox));
      }
      ImageView<PixelMask<Vector2f> > disparity_map_patch = crop(edge_extend(m_disparity_map, ZeroEdgeExtension()),
                                                                 left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      for (int v = 0; v < disparity_map_patch.rows(); ++v)
        for (int u = 0; u < disparity_map_patch.cols(); ++u)
          if ( is_valid(disparity_map_patch(u,v)) )
            remove_mask(disparity_map_patch(u,v)) -= search_range.min();

      switch (m_do_affine_subpixel){

      case 0 : // No Subpixel
        break;
      case 1 : // Parabola Subpixel
        subpixel_correlation_parabola(disparity_map_patch,
                                      left_image_patch,
                                      right_image_patch,
                                      m_kern_width, m_kern_height,
                                      m_do_h_subpixel, m_do_v_subpixel,
                                      m_verbose);
        break;
      case 2: // Bayes EM  Subpixel
        {

          int pyramid_levels = 3;

          // create the pyramid first
          std::vector<ImageView<float> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
          std::vector<BBox2i> regions_of_interest(pyramid_levels);
          left_pyramid[0] = channels_to_planes(left_image_patch);
          right_pyramid[0] = channels_to_planes(right_image_patch);
          regions_of_interest[0] = BBox2i(m_kern_width,m_kern_height,
                                          bbox.width(),bbox.height());

          std::vector<ImageView<PixelMask<Vector2f> > > disparity_map_pyramid(pyramid_levels);
          std::vector<ImageView<PixelMask<Vector2f> > > disparity_map_upsampled(pyramid_levels);
          disparity_map_pyramid[0] = disparity_map_patch;

          // downsample the disparity map and the image pair
          for (int i = 1; i < pyramid_levels; i++) {
            left_pyramid[i] = subsample_img_by_two(left_pyramid[i-1]);
            right_pyramid[i] = subsample_img_by_two(right_pyramid[i-1]);
            disparity_map_pyramid[i] = subsample_disp_map_by_two(disparity_map_pyramid[i-1]);
            regions_of_interest[i] = BBox2i(regions_of_interest[i-1].min()/2, regions_of_interest[i-1].max()/2);
          }

          subpixel_correlation_affine_2d_EM(disparity_map_pyramid[pyramid_levels-1],
                                            left_pyramid[pyramid_levels-1],
                                            right_pyramid[pyramid_levels-1],
                                            m_kern_width, m_kern_height,
                                            regions_of_interest[pyramid_levels-1],
                                            m_do_h_subpixel, m_do_v_subpixel,
                                            m_verbose);
          disparity_map_upsampled[pyramid_levels-1] = disparity_map_pyramid[pyramid_levels-1];

          for (int i = pyramid_levels-2; i >= 0; i--) {

            int up_width = left_pyramid[i].cols();
            int up_height = left_pyramid[i].rows();
            disparity_map_upsampled[i] = upsample_disp_map_by_two(disparity_map_upsampled[i+1], up_width, up_height);

            subpixel_correlation_affine_2d_EM(disparity_map_upsampled[i],
                                              left_pyramid[i],
                                              right_pyramid[i],
                                              m_kern_width, m_kern_height,
                                              regions_of_interest[i],
                                              m_do_h_subpixel, m_do_v_subpixel,
                                              m_verbose);
          }

          disparity_map_patch = disparity_map_upsampled[0];

        }
        break;
      default:
        vw_throw(ArgumentErr() << "Unknown subpixel correlation type: " << m_do_affine_subpixel << ".");
        break;
      }

      // Undo the above adjustment
      for (int v = 0; v < disparity_map_patch.rows(); ++v)
        for (int u = 0; u < disparity_map_patch.cols(); ++u)
          if ( is_valid(disparity_map_patch(u,v)) )
            remove_mask(disparity_map_patch(u,v)) += search_range.min();

      // This may seem confusing, but we must crop here so that the
      // good pixel data is placed into the coordinates specified by
      // the bbox.  This allows rasterize to touch those pixels
      // using the coordinates inside the bbox.  The pixels outside
      // those coordinates are invalid, and they never get accessed.
      return crop(disparity_map_patch, BBox2i(m_kern_width-bbox.min().x(),
                                              m_kern_height-bbox.min().y(),
                                              m_left_image.cols(),
                                              m_left_image.rows() ));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
