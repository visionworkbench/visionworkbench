// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_DISPARITY_MAP_H__
#define __VW_STEREO_DISPARITY_MAP_H__

#include <vw/Core/ProgressCallback.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/Algorithms.h>

// For the PixelDisparity math.
#include <boost/operators.hpp>

namespace vw {

  // Registering the Pixel Disparity type
  template<> struct PixelFormatID<PixelMask<Vector2f> > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<PixelMask<Vector2> > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };

namespace stereo {

  //  get_disparity_range()
  //
  // Determine the range of disparity values present in the disparity map.
  template <class ViewT>
  BBox2 get_disparity_range(ImageViewBase<ViewT> const& disparity_map ) {
    typename UnmaskedPixelType<typename ViewT::pixel_type>::type min, max;
    min_max_pixel_values(disparity_map, min, max);
    return BBox2(min,max);
  }

  //  missing_pixel_image()
  //
  /// Produce a colorized image depicting which pixels in the disparity
  /// map are good pixels, and which are missing (i.e. where no
  /// correlation was found).
  template <class PixelT>
  struct MissingPixelImageFunc: public vw::ReturnFixedType<PixelRGB<uint8> > {
    PixelRGB<uint8> operator() (PixelT const& pix) const {
      if ( is_valid(pix) )
        return PixelRGB<uint8>(200,200,200);
      else
        return PixelRGB<uint8>(255,0,0);
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT, MissingPixelImageFunc<typename ViewT::pixel_type> >
  missing_pixel_image(ImageViewBase<ViewT> &image) {
    return per_pixel_filter(image.impl(), MissingPixelImageFunc<typename ViewT::pixel_type>());
  }

  //  generate_mask()
  //
  /// HAS BEEN REMOVED. See edge_mask() in MaskViews.h

  //  mask_black_pixels()
  //
  /// HAS BEEN REMOVED. See MaskViews.h

  //  disparity_mask()
  //
  //  ......formerly mask()
  //
  /// Given a pair of masks for the left and right images and a
  /// disparity map to be masked, this view will eliminate any pixels
  /// in the disparity map that correspond to locations in the mask
  /// that contain a value of zero.
  template <class PixelT, class MaskViewT>
  struct DisparityMaskFunc: public vw::ReturnFixedType<PixelT>  {

    MaskViewT const& m_left_mask;
    MaskViewT const& m_right_mask;

    DisparityMaskFunc( MaskViewT const& left_mask, MaskViewT const& right_mask) :
      m_left_mask(left_mask), m_right_mask(right_mask) {}

    PixelT operator() (PixelT const& pix, Vector2 const& loc) const {
      if ( !is_valid(pix) ||
           loc[0] < 0 || loc[0] >= m_left_mask.cols() ||
           loc[1] < 0 || loc[1] >= m_left_mask.rows() ||
           loc[0]+pix[0] < 0 || loc[0]+pix[0] >= m_right_mask.cols() ||
           loc[1]+pix[1] < 0 || loc[1]+pix[1] >= m_right_mask.rows() ||
           m_left_mask(vw::int32(loc[0]),vw::int32(loc[1])) == 0 ||
           m_right_mask(vw::int32(loc[0]+pix[0]),vw::int32(loc[1]+pix[1])) == 0 ){
        return PixelT();
      } else
        return pix;
    }
  };

  /// Remove pixels in the disparity map that correspond to locations
  /// where the left or right image mask is invalid. This function
  /// also removes any pixels that fall outside the bounds of
  /// the left or right mask image.
  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT, PixelIndexView, DisparityMaskFunc<typename ViewT::pixel_type,MaskViewT> >
  disparity_mask ( ImageViewBase<ViewT> const& disparity_map,
                   ImageViewBase<MaskViewT> const& left_mask,
                   ImageViewBase<MaskViewT> const& right_mask ) {
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    typedef DisparityMaskFunc<typename ViewT::pixel_type,MaskViewT> func_type;
    typedef BinaryPerPixelView<ViewT, PixelIndexView, func_type > view_type;
    return view_type(disparity_map.impl(),
                     PixelIndexView(disparity_map.impl()),
                     func_type(left_mask.impl(), right_mask.impl()));
  }

  //  disparity_range_mask()
  //
  //  .....formerly remove_invalid_pixels()
  //
  /// Remove pixels from the disparity map that are outside of the
  /// bounds of the original input images.  This happens sometimes
  /// when subpixel interpolation is applied, and this will cause
  /// problems with some camera models (i.e. linear pushbroom) that
  /// are not well defined outside of the bounds of the image.
  template <class PixelT>
  class DisparityRangeMaskFunc: public vw::ReturnFixedType<PixelT>  {

    typedef typename UnmaskedPixelType<PixelT>::type unmasked_type;

    unmasked_type m_min, m_max;

  public:

    DisparityRangeMaskFunc(PixelT const& min, PixelT const& max ) :
    m_min(remove_mask(min)), m_max(remove_mask(max)) {}

    PixelT operator() (PixelT const& pix, Vector2 const& loc) const {
      if ( is_valid(pix) ) {
        if ( loc[0]+pix[0] < m_min[0] || loc[0]+pix[0] >= m_max[0]-1 ||
             loc[1]+pix[1] < m_min[0] || loc[1]+pix[1] >= m_max[1]-1 )
          return PixelT(); // return invalid
      }
      return pix;
    }
  };

  template <class ViewT>
  BinaryPerPixelView<ViewT, PixelIndexView, DisparityRangeMaskFunc<typename ViewT::pixel_type> >
  disparity_range_mask( ImageViewBase<ViewT> &disparity_map,
                        typename ViewT::pixel_type const& min,
                        typename ViewT::pixel_type const& max ) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    typedef DisparityRangeMaskFunc<typename ViewT::pixel_type> func_type;
    typedef BinaryPerPixelView<ViewT, PixelIndexView, func_type> view_type;
    return view_type( disparity_map.impl(),
                      PixelIndexView(disparity_map),
                      func_type( min, max) );
  }


  ///  remove_outliers()
  ///
  /// Remove outliers from a disparity map image using a morphological
  /// erode-like operation.
  template <class PixelT>
  class RemoveOutliersFunc : public ReturnFixedType<PixelT>
  {

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelAccessorView.
    // By maintaining a smart pointer to this small status class, we
    // can change state that is shared between any copies of the
    // RemoveOutliersFunc object and the original.
    struct RemoveOutliersState {
      int rejected_points, total_points;
    };

    int m_half_h_kernel, m_half_v_kernel;
    float m_pixel_threshold;
    float m_rejection_threshold;
    boost::shared_ptr<RemoveOutliersState> m_state;

  public:
    RemoveOutliersFunc(int half_h_kernel, int half_v_kernel, float pixel_threshold, float rejection_threshold) :
      m_half_h_kernel(half_h_kernel), m_half_v_kernel(half_v_kernel),
      m_pixel_threshold(pixel_threshold), m_rejection_threshold(rejection_threshold),
      m_state( new RemoveOutliersState() ) {
      m_state->rejected_points = m_state->total_points = 0;

      VW_ASSERT(half_h_kernel > 0 && half_v_kernel > 0,
                ArgumentErr() << "RemoveOutliersFunc: half kernel sizes must be non-zero.");
    }

    int half_h_kernel() const { return m_half_h_kernel; }
    int half_v_kernel() const { return m_half_v_kernel; }
    float rejection_threshold() const { return m_rejection_threshold; }
    float pixel_threshold() const { return m_pixel_threshold; }
    int rejected_points() const { return m_state->rejected_points; }
    int total_points() const { return m_state->total_points; }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_half_h_kernel, -m_half_v_kernel),
                                             Vector2i(m_half_h_kernel, m_half_v_kernel)); }

    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const {
      m_state->total_points++;

      if (is_valid(*acc)) {
        int matched = 0, total = 0;
        PixelAccessorT row_acc = acc;
        row_acc.advance(-m_half_h_kernel,-m_half_v_kernel);
        for(int yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk) {
          PixelAccessorT col_acc = row_acc;
          for(int xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk) {

            if( is_valid(*col_acc) &&
                fabs((*acc)[0]-(*col_acc)[0]) <= m_pixel_threshold &&
                fabs((*acc)[1]-(*col_acc)[1]) <= m_pixel_threshold) {
              matched++;
            }
            col_acc.next_col();
            total++;
          }
          row_acc.next_row();
        }
        if( ((float)matched/(float)total) < m_rejection_threshold){
          m_state->rejected_points++;
          return typename PixelAccessorT::pixel_type();  //Return invalid pixel
        }
      }
      return *acc;
    }
  };

  // Useful routine for printing how many points have been rejected
  // using a particular RemoveOutliersFunc.
  template <class PixelT>
  inline std::ostream& operator<<(std::ostream& os, RemoveOutliersFunc<PixelT> const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices ("
       << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, RemoveOutliersFunc<typename ViewT::pixel_type> >
  remove_outliers(ImageViewBase<ViewT> const& disparity_map,
                  int half_h_kernel, int half_v_kernel,
                  double pixel_threshold,
                  double rejection_threshold) {
    typedef RemoveOutliersFunc<typename ViewT::pixel_type> func_type;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, func_type > view_type;
    return view_type(edge_extend(disparity_map.impl(), ZeroEdgeExtension()),
                     func_type(half_h_kernel, half_v_kernel,
                               pixel_threshold, rejection_threshold));
  }


  /// Clean up a disparity map.
  ///
  /// You supply the half dimensions of the kernel window.
  ///
  /// Next, you supply the threshold that determines whether a
  /// pixel is considered "close" to its neightbors (in units of
  /// pixels).
  ///
  /// Finally, you supply the percentage of the pixels within the kernel
  /// that must "match" the center pixel if that pixel is to be
  /// considered an inlier. ([0..1.0]).
  template <class ViewT>
  inline UnaryPerPixelAccessorView< EdgeExtensionView< UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>,
                                                                                 RemoveOutliersFunc<typename ViewT::pixel_type> >,
                                                       ZeroEdgeExtension>,
                                    RemoveOutliersFunc<typename ViewT::pixel_type> >
  disparity_clean_up(ImageViewBase<ViewT> const& disparity_map,
                     int h_half_kernel, int v_half_kernel,
                     double pixel_threshold, double rejection_threshold) {
    // Remove outliers first using user specified parameters, and then
    // using a heuristic that isolates single pixel outliers.
    return remove_outliers(remove_outliers(disparity_map.impl(),
                                           h_half_kernel, v_half_kernel,
                                           pixel_threshold, rejection_threshold),
                           1, 1, 1.0, 0.75);
  }


  //  std_dev_image()
  //
  /// Remove pixels from the disparity map that correspond to low
  /// contrast pixels in the original image.
  class StdDevImageFunc : public UnaryReturnTemplateType<PixelTypeFromPixelAccessor>
  {
    int m_kernel_width, m_kernel_height;

  public:
    StdDevImageFunc(int kernel_width, int kernel_height) :
      m_kernel_width(kernel_width), m_kernel_height(kernel_height) {
      VW_ASSERT(m_kernel_width > 0 && m_kernel_height > 0,
                ArgumentErr() << "StdDevImageFunc: kernel sizes must be non-zero.");
    }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_kernel_width/2, -m_kernel_height/2),
                                             Vector2i(m_kernel_width, m_kernel_height)); }

    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const {
      typedef typename PixelAccessorT::pixel_type pixel_type;

      // First pass, compute the mean.
      pixel_type sum = 0;
      PixelAccessorT row_acc = acc;
      row_acc.advance(-m_kernel_width/2,-m_kernel_height/2);
      for(int yk = -m_kernel_height/2; yk <= m_kernel_height/2; ++yk) {
        PixelAccessorT col_acc = row_acc;
        for(int xk = -m_kernel_width/2; xk <= m_kernel_width/2; ++xk) {
          sum += *col_acc;
          col_acc.next_col();
        }
        row_acc.next_row();
      }
      pixel_type mean = sum / (m_kernel_width*m_kernel_height);

      // Second pass, compute the standard deviation using the unbiased
      // estimator.
      sum = 0;
      row_acc = acc;
      row_acc.advance(-m_kernel_width/2,-m_kernel_height/2);
      for(int yk = -m_kernel_height/2; yk <= m_kernel_height/2; ++yk) {
        PixelAccessorT col_acc = row_acc;
        for(int xk = -m_kernel_width/2; xk <= m_kernel_width/2; ++xk) {
          pixel_type diff = *col_acc-mean;
          sum += diff*diff;
          col_acc.next_col();
        }
        row_acc.next_row();
      }
      return sum / (m_kernel_width*m_kernel_height-1);
    }
  };

  template <class ViewT, class EdgeT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, StdDevImageFunc>
  std_dev_image(ImageViewBase<ViewT> const& image,
                int kernel_width, int kernel_height,
                EdgeT edge) {
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, StdDevImageFunc> view_type;
    return view_type(edge_extend(image.impl(), edge),
                     StdDevImageFunc (kernel_width, kernel_height));
  }
  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, StdDevImageFunc>
  std_dev_image(ImageViewBase<ViewT> const& image,
                int kernel_width, int kernel_height) {
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, StdDevImageFunc> view_type;
    return view_type(edge_extend(image.impl(), ZeroEdgeExtension()),
                     StdDevImageFunc (kernel_width, kernel_height));
  }

  // transform_disparities()
  //
  // This Per pixel filter applies an arbitrary transform functor to
  // the pixel coordinates in the secondary image that are encoded by
  // the disparity map.  This is useful for removing the effect of any
  // linear warping for pre-alignment that was performed on the source
  // images prior to correlation.
  template <class TransformT, class PixelT>
  class TransformDisparitiesFunc: public ReturnFixedType<PixelT> {
    TransformT m_trans;

  public:
  TransformDisparitiesFunc(TransformT const& trans) : m_trans(trans) {}

    PixelT operator() (PixelT const& pix, Vector2 const& loc) const {
      Vector2 old_point(loc[0] + pix[0],
                        loc[1] + pix[1]);
      Vector2 new_point = m_trans.reverse(old_point);       // apply the inverse transform
      if ( is_valid(pix) )
        return PixelT(new_point[0] - loc[0],
                      new_point[1] - loc[1]);
      else {
        PixelT result; // invalid
        result[0] = new_point[0]-loc[0];
        result[1] = new_point[1]-loc[1];
        return result;
      }
    }
  };

  template <class ViewT, class TransformT>
  BinaryPerPixelView<ViewT, PixelIndexView, TransformDisparitiesFunc<TransformT, typename ViewT::pixel_type> >
  transform_disparities(ImageViewBase<ViewT> const& disparity_map, TransformT const& transform) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    typedef TransformDisparitiesFunc<TransformT, typename ViewT::pixel_type> func_type;
    typedef BinaryPerPixelView<ViewT, PixelIndexView, func_type > view_type;
    return view_type(disparity_map.impl(),PixelIndexView(disparity_map),
                     func_type(transform));
  }

}}    // namespace vw::stereo


#endif //  _VWDISPARITYMAP_H_

