// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#ifndef __VW_STEREO_DISPARITY_MAP_H__
#define __VW_STEREO_DISPARITY_MAP_H__

#include <vw/Math/Matrix.h>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Core/ProgressCallback.h>

#include <vw/FileIO.h>

// For the PixelDisparity math.
#include <boost/operators.hpp>

namespace vw { 
  
  /// The disparity pixel type for vision workbench has two channels of
  /// disparity (horizontal, vertical) and one alpha channel for recording
  /// which pixels are valid.
  ///
  /// An alpha value of 0 denotes a bad (missing) pixel, and an alpha
  /// channel of 1 denotes a good pixel.
  template <class ChannelT>
  class PixelDisparity : 
    boost::additive< PixelDisparity<ChannelT> >,
    boost::multiplicative< PixelDisparity<ChannelT>, ChannelT >
  {
  private:
    ChannelT m_ch[3];
  public:
    PixelDisparity() { m_ch[0]=m_ch[1]=0; m_ch[2]=1; }
    PixelDisparity( ChannelT const& scalar ) { m_ch[0]=m_ch[1]=scalar; m_ch[2] = 0; }
    PixelDisparity( ChannelT const& h, ChannelT const& v ) { m_ch[0]=h; m_ch[1]=v; m_ch[2]=0; }
  
    template <class OtherT> explicit PixelDisparity( PixelDisparity<OtherT> const& other ) { m_ch[0]=other[0]; m_ch[1]=other[1]; m_ch[2]=other[2];}

    ChannelT& operator[](int i) { return m_ch[i]; }
    ChannelT const& operator[](int i) const { return m_ch[i]; }
    ChannelT& h() { return m_ch[0]; }
    ChannelT const& h() const { return m_ch[0]; }
    ChannelT& v() { return m_ch[1]; }
    ChannelT const& v() const { return m_ch[1]; }
    ChannelT& missing() { return m_ch[2]; }
    ChannelT missing() const { return m_ch[2]; }
    double magnitude() const { return sqrtf((double)m_ch[0]*m_ch[0] + (double)m_ch[1]*m_ch[1]); }
    double magnitude_squared() const { return (double)m_ch[0]*m_ch[0] + (double)m_ch[1]*m_ch[1]; }

    PixelDisparity& operator-() { if (!missing()) {m_ch[0]=-m_ch[0]; m_ch[1]=-m_ch[1];} return *this; }
    PixelDisparity& operator+=( PixelDisparity const& p ) { if (missing() || p.missing()) {m_ch[0]=m_ch[1]=0; m_ch[2] = 1;} 
      else { m_ch[0]+=p[0]; m_ch[1]+=p[1]; } return *this; }
    PixelDisparity& operator-=( PixelDisparity const& p ) { if (missing() || p.missing()) {m_ch[0]=m_ch[1]=0; m_ch[2] = 1;} 
      else { m_ch[0]-=p[0]; m_ch[1]-=p[1]; } return *this; }
    PixelDisparity& operator*=( ChannelT s ) { if (!missing()) { m_ch[0]*=s; m_ch[1]*=s; } return *this; }
    PixelDisparity& operator/=( ChannelT s ) { if (!missing()) { m_ch[0]/=s; m_ch[1]/=s; } return *this; }
    bool operator==(PixelDisparity const& p ) const { 
      return ( missing() == p.missing() && h() == p.h() && v() == p.v() );
    }
  };

  /// \cond INTERNAL
  VW_DECLARE_PIXEL_TYPE(PixelDisparity,3);
  /// \endcond

  template <class ChannelT>
  std::ostream& operator<<( std::ostream& os, PixelDisparity<ChannelT> const& pix ) {
    if (pix.missing()) 
      return os << "Disparity( MISSING )";
    else 
      return os << "Disparity(" << pix.h() << "," << pix.v() << ")";
  }

} // namespace vw

namespace vw {
namespace stereo {
namespace disparity {

  //  get_disparity_range()
  //
  // Determine the range of disparity values present in the disparity map.
  template <class ViewT>
  BBox2 get_disparity_range(ImageViewBase<ViewT> const& disparity_map, int& num_good, bool verbose = false,
                            const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    // Initialize the progress callback
    progress_callback.report_progress(0);

    const ViewT& disparity_map_impl = disparity_map.impl();

    float max_horz_disp = ScalarTypeLimits<float>::lowest();
    float min_horz_disp = ScalarTypeLimits<float>::highest();
    float max_vert_disp = ScalarTypeLimits<float>::lowest();
    float min_vert_disp = ScalarTypeLimits<float>::highest();

    // Find the max/min disparity values
    num_good = 0;
    for (int32 j = 0; j < disparity_map_impl.rows(); j++) {

      // Update the progress callback.
      if (progress_callback.abort_requested()) 
        vw_throw( Aborted() << "Aborted by ProgressCallback" );
      progress_callback.report_progress((float)j/disparity_map_impl.rows());

      for (int32 i = 0; i < disparity_map_impl.cols(); i++) {
        if ( !disparity_map_impl(i,j).missing() ) {
          max_horz_disp = disparity_map_impl(i,j).h() > max_horz_disp ? disparity_map_impl(i,j).h() : max_horz_disp;
          min_horz_disp = disparity_map_impl(i,j).h() < min_horz_disp ? disparity_map_impl(i,j).h() : min_horz_disp;
          max_vert_disp = disparity_map_impl(i,j).v() > max_vert_disp ? disparity_map_impl(i,j).v() : max_vert_disp;
          min_vert_disp = disparity_map_impl(i,j).v() < min_vert_disp ? disparity_map_impl(i,j).v() : min_vert_disp;
          num_good++;
        }
      }
    }
    progress_callback.report_finished();

    if (num_good == 0) {
      if (verbose)
        vw_out(WarningMessage, "stereo") << "Disparity range -- disparity map had zero good pixels.";
      return BBox2(0,0,0,0);
    }
    
    if (verbose) 
      vw_out(InfoMessage, "stereo") << "Disparity range -- Horizontal: [" << min_horz_disp << ", " << max_horz_disp 
                                    << "]   Vertical: [" << min_vert_disp << ", " << max_vert_disp << "]  ("<< num_good << " good)\n"; 
    return BBox2(Vector2(min_horz_disp, min_vert_disp),Vector2(max_horz_disp, max_vert_disp));
  }


  template <class ViewT>
  BBox2 get_disparity_range(ImageViewBase<ViewT> const& disparity_map, bool verbose = false) {
    int num_good = 0;
    return get_disparity_range(disparity_map.impl(), num_good, verbose);
  }

  //  missing_pixel_image()
  //
  /// Produce a colorized image depicting which pixels in the disparity
  /// map are good pixels, and which are missing (i.e. where no
  /// correlation was found).
  struct MissingPixelImageFunc: public vw::ReturnFixedType<PixelRGB<uint8> > {
    PixelRGB<uint8> operator() (PixelDisparity<float> const& pix) const {
      if ( !pix.missing() ) 
        return PixelRGB<uint8>(200,200,200);
      else
        return PixelRGB<uint8>(255,0,0);
    }
  };
    
  template <class ViewT>
  UnaryPerPixelView<ViewT, MissingPixelImageFunc> 
  missing_pixel_image(ImageViewBase<ViewT> &disparity_map) {
    return per_pixel_filter(disparity_map.impl(), MissingPixelImageFunc());
  }


  //  generate_mask()
  // 
  /// Masks all of the black pixels along the edges of an image.  This
  /// algorithm "eats away" at the pixels on all four sides of the
  /// image; masking pixels until it encounters a non-black pixel.
  /// 
  /// You can supply an optional buffer argument that will mask some
  /// of the good data bordering the image as well.  The width of this
  /// border is set using the additional_border_width argument (in
  /// units of pixels).
  template <class ViewT>
  vw::ImageView<bool> generate_mask(vw::ImageViewBase<ViewT> const& input_image,
                                    int32 additional_border_width = 0) {

    typedef typename ViewT::pixel_type pixel_type;
    vw::ImageView<bool> mask(input_image.impl().cols(), input_image.impl().rows(), 1);
    int32 i, j;

    vw::fill(mask, true);

    for (i = 0; i < input_image.impl().cols(); i++) {
      // Search from the left side of the image for black pixels 
      j = 0;
      while ( j < input_image.impl().rows() && input_image.impl()(i,j) == pixel_type() ) {
        mask(i,j) = false;
        j++;
      }

      // Mask additional "buffer" pixels
      int j_start = j;
      while ( j < input_image.impl().rows() && j - j_start < additional_border_width ) {
        mask(i,j) = false;
        j++;
      }
    
      // Search from the right side of the image for black pixels 
      j = input_image.impl().rows() - 1;
      while ( j > 0 && input_image.impl()(i,j) == pixel_type() ) {
        mask(i,j) = false;
        j--;
      }

      // Mask additional "buffer" pixels
      int j_end = j;
      while ( j > 0 && j_end - j < additional_border_width) {
        mask(i,j) = false;
        j--;
      }
    }

    for (j = 0; j < input_image.impl().rows(); j++) {
      // Search from the top side of the image for black pixels 
      i = 0;
      while ( i < input_image.impl().cols() && input_image.impl()(i,j) == pixel_type() ) {
        mask(i,j) = false;
        i++;
      }

      // Mask additional "buffer" pixels
      int i_start = i;
      while ( i < input_image.impl().cols() && i - i_start < additional_border_width) {
        mask(i,j) = false;
        i++;
      }
    
      // Search from the bottom side of the image for black pixels 
      i = input_image.impl().cols() - 1;
      while ( i > 0 && input_image.impl()(i,j) == pixel_type()) {
        mask(i,j) = false;
        i--;
      }

      // Mask additional "buffer" pixels
      int i_end = i;
      while ( i > 0 && i_end - i < additional_border_width) {
        mask(i,j) = false;
        i--;
      }
    }

    return mask;
  }

  //  mask_black_pixels()
  template <class ViewT>
  void mask_black_pixels(vw::ImageViewBase<ViewT> const& input_image,
                         vw::ImageView<uint8> &mask_image) {

    typedef typename ViewT::pixel_type pixel_type;
    
    VW_ASSERT(input_image.impl().cols() == mask_image.cols() && input_image.impl().rows() == mask_image.rows(),
              ArgumentErr() << "mask_black_pixels() : Input image and mask images do not have the same dimensions.");

    for (int j=0; j < mask_image.rows(); ++j) {
      for (int i=0; i < mask_image.cols(); ++i) {
        if (input_image.impl()(i,j) == pixel_type()) 
          mask_image(i,j) = 0;
      }
    }
  }

  //  mask()
  //
  /// Given a pair of masks for the left and right images and a
  /// disparity map to be masked, this view will eliminate any pixels
  /// in the disparity map that correspond to locations in the mask
  /// that contain a value of zero.
  template <class MaskViewT>
  struct DisparityMaskFunc: public vw::ReturnFixedType<PixelDisparity<float> >  {

    MaskViewT m_left_mask; 
    MaskViewT m_right_mask;

    DisparityMaskFunc( MaskViewT const& left_mask, MaskViewT const& right_mask) :
      m_left_mask(left_mask), m_right_mask(right_mask) {}

    PixelDisparity<float> operator() (PixelDisparity<float> const& pix, Vector3 const& loc) const {
      if ( pix.missing() ||                                               // If already a missing pixel
           loc[0] < 0 || loc[0] >= m_left_mask.cols() ||                  //
           loc[1] < 0 || loc[1] >= m_left_mask.rows() ||                  // or outside of bounds of 
           loc[0]+pix.h() < 0 || loc[0]+pix.h() >= m_right_mask.cols() || // the left or right image
           loc[1]+pix.v() < 0 || loc[1]+pix.v() >= m_right_mask.rows() || //
           !(m_left_mask(int(loc[0]), int(loc[1]))) ||                    // or the pixel is masked
           !(m_right_mask(int(loc[0]+pix.h()), int(loc[1]+pix.v()))) ) {
        return PixelDisparity<float>();                                   // then set to missing pixel value
      } else {
        return pix;
      } 

    }
  };

  /// Remove pixels in the disparity map that correspond to locations
  /// where the left or right image mask contains zeros.  This
  /// function also removes any pixels that fall outside the bounds of
  /// the left or right mask image.
  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT, PixelIndex3View, DisparityMaskFunc<MaskViewT> > mask(ImageViewBase<ViewT> const& disparity_map, 
                                                                                 ImageViewBase<MaskViewT> const& left_mask, 
                                                                                 ImageViewBase<MaskViewT> const& right_mask) {
    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    return BinaryPerPixelView<ViewT, PixelIndex3View, DisparityMaskFunc<MaskViewT> >(disparity_map.impl(), 
                                                                                     PixelIndex3View(disparity_map),
                                                                                     DisparityMaskFunc<MaskViewT>(left_mask.impl(), right_mask.impl()));
  }



  //  remove_invalid_pixels()
  //
  /// Remove pixels from the disparity map that are outside of the
  /// bounds of the original input images.  This happens sometimes
  /// when subpixel interpolation is applied, and this will cause
  /// problems with some camera models (i.e. linear pushbroom) that
  /// are not well defined outside of the bounds of the image.
  struct InvalidPixelsFunc: public vw::ReturnFixedType<PixelDisparity<float> >  {

    int m_width, m_height;

    InvalidPixelsFunc(int right_image_width, int right_image_height) : 
      m_width(right_image_width), m_height(right_image_height) {}

    PixelDisparity<float> operator() (PixelDisparity<float> const& pix, Vector3 const& loc) const {
      if ( !pix.missing() ) {
        if ( loc[0]+pix.h() < 0 || loc[0]+pix.h() >= m_width-1 || 
             loc[1]+pix.v() < 0 || loc[1]+pix.v() >= m_height-1) {
          return PixelDisparity<float>(); // Set to missing pixel
        } 
      }
      return pix;
    }
  };
    
  template <class ViewT>
  BinaryPerPixelView<ViewT, PixelIndex3View, InvalidPixelsFunc> 
  remove_invalid_pixels(ImageViewBase<ViewT> &disparity_map, int right_image_width, int right_image_height) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    return BinaryPerPixelView<ViewT, PixelIndex3View, InvalidPixelsFunc>(disparity_map.impl(), 
                                                                        PixelIndex3View(disparity_map),
                                                                        InvalidPixelsFunc(right_image_width,
                                                                                          right_image_height) );
  }
  



  //  remove_border_pixels()
  //
  /// Remove pixels from the disparity map that lie along the borders.
  struct BorderPixelsFunc: public vw::ReturnFixedType<PixelDisparity<float> >  {

    int m_border_size, m_image_width, m_image_height;

    BorderPixelsFunc(int border_size, int image_width, int image_height) : 
      m_border_size(border_size), m_image_width(image_width), m_image_height(image_height) {}

    PixelDisparity<float> operator() (PixelDisparity<float> const& pix, Vector3 const& loc) const {
      if ( loc[0] < m_border_size || loc[1] < m_border_size ||
           loc[0] > m_image_width-m_border_size ||
           loc[1] > m_image_height-m_border_size) {
        return PixelDisparity<float>(); // Set to missing pixel
      } 
      return pix;
    }
  };
    
  template <class ViewT>
  BinaryPerPixelView<ViewT, PixelIndex3View, BorderPixelsFunc> 
  remove_border_pixels(ImageViewBase<ViewT> const& disparity_map, int border_size) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    return BinaryPerPixelView<ViewT, PixelIndex3View, BorderPixelsFunc>(disparity_map.impl(), 
                                                                       PixelIndex3View(disparity_map),
                                                                       BorderPixelsFunc(border_size,
                                                                                         disparity_map.cols(),
                                                                                         disparity_map.rows()) );
  }
  
  //  remove_outliers()
  // 
  /// Remove outliers from a disparity map image using a morpholical
  /// approach. 
  class RemoveOutliersFunc : public ReturnFixedType<PixelDisparity<float> > 
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

      if (!(*acc).missing()) {        
        int matched = 0, total = 0; 
        PixelAccessorT row_acc = acc;
        row_acc.advance(-m_half_h_kernel,-m_half_v_kernel); 
        for(int yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk) {
          PixelAccessorT col_acc = row_acc;
          for(int xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk) {

            if( !(*col_acc).missing() &&
                fabs((*acc).h()-(*col_acc).h()) <= m_pixel_threshold &&
                fabs((*acc).v()-(*col_acc).v()) <= m_pixel_threshold) {
              matched++;
            }
            col_acc.next_col();
            total++;
          }
          row_acc.next_row();
        }
        if( ((float)matched/(float)total) < m_rejection_threshold){
          m_state->rejected_points++;
          return typename PixelAccessorT::pixel_type();
        }
      } 
      return *acc;
    }
  };

  // Useful routine for printing how many points have been rejected
  // using a particular RemoveOutliersFunc.
  inline std::ostream& operator<<(std::ostream& os, RemoveOutliersFunc const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices (" << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, RemoveOutliersFunc> remove_outliers(ImageViewBase<ViewT> const& disparity_map,
                                                                                                      int half_h_kernel, int half_v_kernel,
                                                                                                      double pixel_threshold,
                                                                                                      double rejection_threshold) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, RemoveOutliersFunc>(edge_extend(disparity_map.impl(), ZeroEdgeExtension()),
                                                                                               RemoveOutliersFunc (half_h_kernel, 
                                                                                                                   half_v_kernel, 
                                                                                                                   pixel_threshold, 
                                                                                                                   rejection_threshold));
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
  inline UnaryPerPixelAccessorView<EdgeExtensionView<UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, 
                                                                               RemoveOutliersFunc>, 
                                                     ZeroEdgeExtension>, RemoveOutliersFunc>  
  clean_up(ImageViewBase<ViewT> const& disparity_map,
           int h_half_kernel, int v_half_kernel, 
           double pixel_threshold, double rejection_threshold) {
    
    // Remove outliers first using user specified parameters, and then
    // using a heuristic that isolates single pixel outliers.
    return remove_outliers(remove_outliers(disparity_map.impl(), 
                                           h_half_kernel, v_half_kernel, 
                                           pixel_threshold, rejection_threshold),
                           1, 1, 1.0, 0.75);
  }





  //  low_contrast_filter()
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
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, StdDevImageFunc> std_dev_image(ImageViewBase<ViewT> const& image,
                                                                                                          int kernel_width, int kernel_height, 
                                                                                                          EdgeT edge) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, StdDevImageFunc>(edge_extend(image.impl(), edge),
                                                                                                     StdDevImageFunc (kernel_width, kernel_height));
  }
  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, StdDevImageFunc> std_dev_image(ImageViewBase<ViewT> const& image,
                                                                                                          int kernel_width, int kernel_height) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, StdDevImageFunc>(edge_extend(image.impl(), ZeroEdgeExtension()),
                                                                                                     StdDevImageFunc (kernel_width, kernel_height));
  }

  class LessThanThresholdFunc: public vw::ReturnFixedType<bool> {
    double m_threshold;
  public:
    LessThanThresholdFunc(double threshold) : m_threshold(threshold) {}

    template <class PixelT>
    bool operator() (PixelT const& pix) const {
      return pix < m_threshold;
    }
  };
    
  template <class ViewT> UnaryPerPixelView<ViewT, LessThanThresholdFunc> 
  less_than_threshold(ImageViewBase<ViewT> const& image, double threshold) {
    return per_pixel_filter(image.impl(), LessThanThresholdFunc(threshold));
  }
  // transform_disparities()
  //
  // This Per pixel filter applies an arbitrary transform functor to
  // the pixel coordinates in the secondary image that are encoded by
  // the disparity map.  This is useful for removing the effect of any
  // linear warping for pre-alignment that was performed on the source
  // images prior to correlation.
  template <class TransformT>
  class TransformDisparitiesFunc: public ReturnFixedType<PixelDisparity<float> > {
    TransformT m_trans;

  public:
    TransformDisparitiesFunc(TransformT const& trans) : m_trans(trans) {}
    
    PixelDisparity<float> operator() (PixelDisparity<float> const& pix, Vector3 const& loc) const {
      if ( !pix.missing() ) {
        
        Vector2 old_point(loc[0] + pix.h(),
                          loc[1] + pix.v());
        Vector2 new_point = m_trans.reverse(old_point);       // apply the inverse transform
        return PixelDisparity<float>(new_point[0] - loc[0],
                                     new_point[1] - loc[1]);    
      } else {
        return pix;
      }
    }
  };
  
  template <class ViewT, class TransformT>
  BinaryPerPixelView<ViewT, PixelIndex3View, TransformDisparitiesFunc<TransformT> > 
  transform_disparities(ImageViewBase<ViewT> const& disparity_map, TransformT const& transform) {
    
    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    return BinaryPerPixelView<ViewT, PixelIndex3View, TransformDisparitiesFunc<TransformT> >(disparity_map.impl(), 
                                                                                             PixelIndex3View(disparity_map),
                                                                                             TransformDisparitiesFunc<TransformT>(transform));
  }

  

//   template <class ViewT>
//   void sparse_disparity_filter(ImageViewBase<ViewT> const& disparity_map, 
//                                float blur_stddev, float rejection_threshold) {
//     ImageViewRef<PixelGray<float> > test_image = disparity::missing_pixel_image(disparity_map);
//     //    write_image("test-a.png", test_image);
//     ImageView<float> blurred_image = gaussian_filter(select_channel(test_image,0), blur_stddev);
//     //    write_image("test-b.png", blurred_image);
//     ImageView<float> threshold_image = threshold(blurred_image, rejection_threshold);
//     //    write_image("test-c.png", threshold_image);
    
//     for (int i = 0; i < disparity_map.cols(); i++) 
//       for (int j = 0; j < disparity_map.rows(); j++) 
//         if (threshold_image(i,j) == 0)
//           disparity_map(i,j) = PixelDisparity<ChannelT>(); // Set to missing pixel
    
//     std::cout << " done.\n";
//   }  
  
} // namespace disparity
  
}}    // namespace vw::stereo


#endif //  _VWDISPARITYMAP_H_

