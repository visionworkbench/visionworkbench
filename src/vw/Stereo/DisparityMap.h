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


#ifndef __VW_STEREO_DISPARITY_MAP_H__
#define __VW_STEREO_DISPARITY_MAP_H__

#include <vw/vw_config.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Statistics.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/Transform.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>

#include <ostream>

// For the PixelDisparity math.
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/type_traits/remove_reference.hpp>

namespace vw {

namespace stereo {

  //  get_disparity_range()
  //
  /// Determine the range of disparity values present in the disparity map.
  template <class ViewT>
  BBox2f get_disparity_range(ImageViewBase<ViewT> const& disparity_map ) {
    typedef typename UnmaskedPixelType<typename ViewT::pixel_type>::type accum_type;
    PixelAccumulator<EWMinMaxAccumulator<accum_type>> accumulator;

    // TODO(oalexan1): It looks as if this accumulator does not skip
    // invalid pixels.
    for_each_pixel(disparity_map, accumulator);

    // This only checks if any pixels were passed in, not if they were
    // valid or not.
    if (!accumulator.is_valid())
      return BBox2f(0,0,0,0);
    
    return BBox2f(accumulator.minimum(), accumulator.maximum());
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
  missing_pixel_image(ImageViewBase<ViewT> const& image) {
    return per_pixel_filter(image.impl(), MissingPixelImageFunc<typename ViewT::pixel_type>());
  }

  ///  disparity_mask()
  ///
  ///  ......formerly mask()
  ///
  /// Given a pair of masks for the left and right images and a
  /// disparity map to be masked, this view will eliminate any pixels
  /// in the disparity map that correspond to locations in the mask
  /// that contain a value of zero.
  template <class ViewT, class MaskView1T, class MaskView2T>
  class DisparityMaskView : public ImageViewBase<DisparityMaskView<ViewT, MaskView1T, MaskView2T> > {
    ViewT      m_input_view;
    MaskView1T m_mask1_view;
    MaskView2T m_mask2_view;
    BBox2i     m_search_range;

    template <class MaskPT>
    bool has_valid_mask( ImageView<MaskPT> const& mask ) const {
      typedef typename ImageView<MaskPT>::pixel_accessor Pacc;
      Pacc row_acc =
        mask.origin();
      for ( int32 r = mask.rows(); r; --r ) {
        Pacc col_acc = row_acc;
        for ( int32 c = mask.cols(); c; --c ) {
          if ( *col_acc )
            return true;
          col_acc.next_col();
        }
        row_acc.next_row();
      }
      return false;
    }

  public:
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename ViewT::pixel_type result_type;
    typedef ProceduralPixelAccessor<DisparityMaskView> pixel_accessor;

    DisparityMaskView( ImageViewBase<ViewT>      const& image,
                       ImageViewBase<MaskView1T> const& mask1,
                       ImageViewBase<MaskView2T> const& mask2,
                       BBox2i const& search_range = BBox2i() ) :
      m_input_view(image.impl()), m_mask1_view(mask1.impl()), m_mask2_view(mask2.impl()), m_search_range( search_range ) {
      VW_DEBUG_ASSERT( image.impl().cols() == mask1.impl().cols() &&
                       image.impl().rows() == mask1.impl().rows(),
                       ArgumentErr() << "disparity_mask: input and left mask are not same dimensions." );
    }

    // Standard required ImageView interface
    inline int32 cols  () const { return m_input_view.cols  (); }
    inline int32 rows  () const { return m_input_view.rows  (); }
    inline int32 planes() const { return m_input_view.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline result_type operator()( int32 i, int32 j, int32 p = 0 ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      if ( i < 0 || j < 0 || p < 0 ||
           i > m_mask1_view.cols() || j > m_mask1_view.rows() ||
           p > m_mask1_view.planes() ) {
        vw_throw( ArgumentErr() << "DisparityMaskView::operator() - invalid index [" << i << " "
                  << j << " " << p << "]" );
      }
#endif
      if ( !m_mask1_view(i,j,p) )
        return result_type();
      result_type disparity =
        m_input_view(i,j,p);
      if ( !is_valid(disparity) )
        return result_type();
      if ( i+disparity[0] < 0 || i+disparity[0] >= m_mask2_view.cols() ||
           j+disparity[1] < 0 || j+disparity[1] >= m_mask2_view.rows() ||
           m_mask2_view(i+disparity[0],j+disparity[1],p) == 0 )
        return result_type();
      return disparity;
    }

    /// Block rasterization section that does actual work
    typedef DisparityMaskView<CropView<ImageView<typename     ViewT::pixel_type > >, 
                              CropView<ImageView<typename MaskView1T::pixel_type> >, 
                              CropView<ImageView<typename MaskView2T::pixel_type> > > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {
      typedef typename MaskView1T::pixel_type Mask1PT;
      typedef typename MaskView2T::pixel_type Mask2PT;
      typedef typename ViewT::pixel_type      ViewPT;

      // Prerasterize Mask1 as we'll always be using it. I'm not using
      // the prerasterize type as I want to make sure I can write to
      // this memory space.
      CropView<ImageView<Mask1PT> > mask1_preraster =
        crop(ImageView<Mask1PT>(crop(m_mask1_view,bbox)),-bbox.min().x(),-bbox.min().y(),cols(),rows());

      // Check to see if mask1 is even occupied, if not let's not
      // prerasterize anything and let this view just return a blank
      // disparity map.
      if ( !has_valid_mask(mask1_preraster.child()) )
        return prerasterize_type( crop(ImageView<ViewPT>(0,0),0,0,
                                       cols(),rows()),
                                  mask1_preraster,
                                  crop(ImageView<Mask2PT>(0,0),0,0,0,0) );

      // We have bbox so now we know what sections of the
      // right mask to prerasterize.
      if ( m_search_range != BBox2i() ) {
        BBox2i right_bbox = bbox;
        right_bbox.min() += m_search_range.min();
        right_bbox.max() += m_search_range.max();
        right_bbox.crop( bounding_box(m_mask2_view) );

        CropView<ImageView<Mask2PT> > mask2_preraster =
          crop( ImageView<Mask2PT>( crop(m_mask2_view,right_bbox) ),
                -right_bbox.min().x(), -right_bbox.min().y(),
                m_mask2_view.cols(), m_mask2_view.rows() );
        if ( !has_valid_mask(mask2_preraster.child()) ) {
          // It appears the right mask is completely empty. In order
          // to make this view early exit, we'll set the left mask to entire zero.
          fill( mask1_preraster.child(), Mask1PT() );
          return prerasterize_type( crop(ImageView<ViewPT>(0,0),0,0,cols(),rows()),
                                    mask1_preraster,
                                    crop(ImageView<Mask2PT>(0,0),0,0,0,0) );
        }

        // Actually rasterize the disparity
        return prerasterize_type( crop(ImageView<ViewPT>(crop(m_input_view,bbox)),
                                       -bbox.min().x(), -bbox.min().y(), cols(), rows() ),
                                  mask1_preraster, mask2_preraster );
      }

      // If we don't know the search range, prerasterize the
      // disparity and work it out by reading the data.
      CropView<ImageView<ViewPT> > disp_preraster =
        crop( ImageView<ViewPT>( crop(m_input_view, bbox) ),
              -bbox.min().x(), -bbox.min().y(), cols(), rows() );
      typedef typename UnmaskedPixelType<ViewPT>::type accum_t;
      PixelAccumulator<EWMinMaxAccumulator<accum_t> > accumulator;
      for_each_pixel( disp_preraster.child(), accumulator );
      if ( !accumulator.is_valid() )
        return prerasterize_type( disp_preraster,
                                  mask1_preraster,
                                  crop(ImageView<Mask2PT>(0,0),0,0,0,0) );
      accum_t input_min = accumulator.minimum();
      accum_t input_max = accumulator.maximum();
      BBox2i preraster(bbox.min() + floor(Vector2f(input_min[0],input_min[1])),
                       bbox.max() + ceil(Vector2f(input_max[0],input_max[1])) );
      preraster.crop( bounding_box( m_mask2_view ) );
      // The bottom might be a little confusing. We're rasterizing the
      // part of mask2 that we know the input disparity will actually touch.
      return prerasterize_type( disp_preraster,
                                mask1_preraster,
                                crop( ImageView<Mask2PT>( crop(m_mask2_view,preraster) ),
                                      -preraster.min().x(), -preraster.min().y(),
                                      m_mask2_view.cols(),
                                      m_mask2_view.rows() ) );
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize( prerasterize(bbox), dest, bbox ); }
  };

  template <class ViewT, class MaskView1T, class MaskView2T>
  DisparityMaskView<ViewT, MaskView1T, MaskView2T>
  disparity_mask( ImageViewBase<ViewT>      const& disparity_map,
                  ImageViewBase<MaskView1T> const& left_mask,
                  ImageViewBase<MaskView2T> const& right_mask,
                  BBox2i const& search_range = BBox2i() ) {
    return DisparityMaskView<ViewT,MaskView1T,MaskView2T>(disparity_map.impl(),left_mask.impl(),right_mask.impl(), search_range);
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
  BinaryPerPixelView<ViewT, PerPixelIndexView<VectorIndexFunctor>, DisparityRangeMaskFunc<typename ViewT::pixel_type> >
  disparity_range_mask( ImageViewBase<ViewT> const& disparity_map,
                        typename ViewT::pixel_type const& min,
                        typename ViewT::pixel_type const& max ) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    typedef DisparityRangeMaskFunc<typename ViewT::pixel_type> func_type;
    typedef BinaryPerPixelView<ViewT, PerPixelIndexView<VectorIndexFunctor>, func_type> view_type;
    return view_type( disparity_map.impl(),
                      pixel_index_view(disparity_map),
                      func_type( min, max) );
  }

  // ================================================================================
  // Start of outlier removal functions.


  /// Remove outliers from a disparity map image.

  // Method 1: Use some thresholds and counting.
  /// You supply the half dimensions of the kernel window.
  ///
  /// Next, you supply the threshold that determines whether a
  /// pixel is considered "close" to its neighbors (in units of pixels).
  ///
  /// Finally, you supply the percentage of the pixels within the kernel
  /// that must "match" the center pixel if that pixel is to be
  /// considered an inlier. ([0..1.0]).

  template <class PixelT>
  class RmOutliersUsingThreshFunc : public ReturnFixedType<PixelT>
  {

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelAccessorView.
    // By maintaining a smart pointer to this small status class, we
    // can change state that is shared between any copies of the
    // RmOutliersUsingThreshFunc object and the original.
    struct RmOutliersState {
      int32 rejected_points, total_points;
    };

    int32 m_half_h_kernel, m_half_v_kernel;
    double m_pixel_threshold;
    double m_rejection_threshold;
    boost::shared_ptr<RmOutliersState> m_state;

  public:
    RmOutliersUsingThreshFunc(int32 half_h_kernel, int32 half_v_kernel, double pixel_threshold, double rejection_threshold) :
      m_half_h_kernel(half_h_kernel), m_half_v_kernel(half_v_kernel),
      m_pixel_threshold(pixel_threshold), m_rejection_threshold(rejection_threshold),
      m_state( new RmOutliersState() ) {
      m_state->rejected_points = m_state->total_points = 0;

      VW_ASSERT(half_h_kernel > 0 && half_v_kernel > 0,
                ArgumentErr() << "RmOutliersFunc: half kernel sizes must be non-zero.");
    }

    int32  half_h_kernel      () const { return m_half_h_kernel;          }
    int32  half_v_kernel      () const { return m_half_v_kernel;          }
    double rejection_threshold() const { return m_rejection_threshold;    }
    double pixel_threshold    () const { return m_pixel_threshold;        }
    int32  rejected_points    () const { return m_state->rejected_points; }
    int32  total_points       () const { return m_state->total_points;    }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_half_h_kernel, -m_half_v_kernel),
                                             Vector2i( m_half_h_kernel,  m_half_v_kernel)); }

    /// Returns either a copy of the input pixel (PASS) or an empty invalid pixel (FAIL)
    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const {
      m_state->total_points++;

      if (is_valid(*acc)) {
        int32 matched = 0, total = 0;
        PixelAccessorT row_acc = acc;
        row_acc.advance(-m_half_h_kernel,-m_half_v_kernel);
        for(int32 yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk) {
          PixelAccessorT col_acc = row_acc;
          for(int32 xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk) {

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
        if( ((double)matched/(double)total) < m_rejection_threshold){
          m_state->rejected_points++;
          return typename PixelAccessorT::pixel_type();  //Return invalid pixel
        }
      }
      return *acc;
    }
  };

  // Useful routine for printing how many points have been rejected
  // using a particular RmOutliersFunc.
  template <class PixelT>
  inline std::ostream&
  operator<<(std::ostream& os, RmOutliersUsingThreshFunc<PixelT> const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices ("
       << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  /// Wrapper to apply the RmOutliersUsingThreshFunc functor to an entire image.
  /// - Using this call takes care of inserting the required EdgeExtension view
  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, 
                            RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >
  rm_outliers_using_thresh(ImageViewBase<ViewT> const& disparity_map,
                           int32 half_h_kernel, int32 half_v_kernel,
                           double pixel_threshold,
                           double rejection_threshold) {
    typedef RmOutliersUsingThreshFunc<typename ViewT::pixel_type> func_type;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, func_type > view_type;
    return view_type(edge_extend(disparity_map.impl(), ConstantEdgeExtension()),
                     func_type(half_h_kernel, half_v_kernel,
                               pixel_threshold, rejection_threshold));
  }

  /// Remove speckle outliers first using user specified parameters, and then
  /// using a heuristic that isolates single pixel outliers.
  /// - Rejection_threshold is percentage of neighbors that must be within pixel_threshold distance
  ///   of the candidate pixel.
  /// - This function is essentially a two-pass version of rm_outliers_using_thresh with one set of
  ///   parameters being hard-coded.
  template <class ViewT>
  inline UnaryPerPixelAccessorView< UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>,
                                                              RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >, 
                                                              RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >
  disparity_cleanup_using_thresh(ImageViewBase<ViewT> const& disparity_map,
                                 int32 h_half_kernel, int32 v_half_kernel,
                                 double pixel_threshold,
                                 double rejection_threshold) {
    typedef RmOutliersUsingThreshFunc<typename ViewT::pixel_type> func_type; // Functor definition
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, func_type > inner_type;
    typedef UnaryPerPixelAccessorView<inner_type, func_type > outer_type; // Second application of func_type
    return outer_type(rm_outliers_using_thresh
                       (disparity_map.impl(),
                        h_half_kernel, v_half_kernel,
                        pixel_threshold, rejection_threshold
                       ),
                      //func_type_thresh( 1, 1, 1.0, 0.75 ) // overly restrictive
                      // at least 1 neighbor must be within 3 pixels
                      func_type( 1, 1, 3.0, 0.20 )
                      );
  }

  // Method 2: Compare the central point to the mean of the neighbors (Ara's method).
  template <class PixelT>
  class RmOutliersUsingMeanFunc : public ReturnFixedType<PixelT>{

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelAccessorView.
    // By maintaining a smart pointer to this small status class, we can
    // change state that is shared between any copies of the
    // RmOutliersFunc object and the original.
    struct RmOutliersState {
      int32 rejected_points, total_points;
    };

    int32 m_half_h_kernel, m_half_v_kernel;
    double m_max_mean_diffSq;
    boost::shared_ptr<RmOutliersState> m_state;

  public:
  
    RmOutliersUsingMeanFunc(int32 half_h_kernel, int32 half_v_kernel,
                            double max_mean_diff) :
      m_half_h_kernel(half_h_kernel), m_half_v_kernel(half_v_kernel),
      m_max_mean_diffSq(max_mean_diff*max_mean_diff),
      m_state( new RmOutliersState() ) {
      m_state->rejected_points = m_state->total_points = 0;
    
      VW_ASSERT(half_h_kernel > 0 && half_v_kernel > 0,
                ArgumentErr() << "RmOutliersUsingMeanFunc: half kernel sizes must be non-zero.");
    }

    int32 half_h_kernel() const { return m_half_h_kernel; }
    int32 half_v_kernel() const { return m_half_v_kernel; }
    double max_mean_diff() const { return sqrt(m_max_mean_diffSq); }
    int32 rejected_points() const { return m_state->rejected_points; }
    int32 total_points() const { return m_state->total_points; }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_half_h_kernel, -m_half_v_kernel),
                                             Vector2i(m_half_h_kernel, m_half_v_kernel)); }

    // This function will be called for each pixel in the image to be filtered
    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const{
      m_state->total_points++; // Accumulate number of points operated on
    
      // Quit immediately if the current point is already invalid
      if (!is_valid(*acc))
        return *acc;

      // Remove gross outliers: Find the 75th percentile largest disparity
      // magnitude, and multiply it by 2. Any disparity with magnitude
      // larger than this will be thrown out. 
      // E.g., if magnitudes are 1, 2, 3, 4, 5, 6, 7, 8, 1e10,
      // 2 * 75%th magnitude is 2*7, so 1e10 is thrown out.
      PixelAccessorT row_acc0 = acc;
      std::vector<double> len;
      row_acc0.advance(-m_half_h_kernel,-m_half_v_kernel); // Move to top left kernel
      for(int32 yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk){
        PixelAccessorT col_acc = row_acc0; // Start column iterator at left of current row
        for(int32 xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk){
          if(is_valid(*col_acc)){
            double val = std::abs((*col_acc)[0]) + std::abs((*col_acc)[1]);
            len.push_back(val);
          }
          col_acc.next_col(); // Advance to next column
        }
        row_acc0.next_row(); // Advance to next row
      }
      std::sort(len.begin(), len.end());
      double cutoff = 0.0;
      if (!len.empty()) cutoff = 2.0*len[(int)(0.75*len.size())];

      // Compute the mean
      size_t matched = 0, total = 0;
      double meanX=0, meanY=0;
      PixelAccessorT row_acc = acc; // position at the center
      row_acc.advance(-m_half_h_kernel,-m_half_v_kernel); // Move to top left kernel
      for(int32 yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk){
        PixelAccessorT col_acc = row_acc; // Start column iterator at left of current row
        for(int32 xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk){
          if(is_valid(*col_acc)){
            double val = std::abs((*col_acc)[0]) + std::abs((*col_acc)[1]);
            if (val > cutoff) continue; // unreasonably large disparity
            meanX += (*col_acc)[0];
            meanY += (*col_acc)[1];
            matched++; // Total number of valid pixels
          }
          col_acc.next_col(); // Advance to next column
          total++; // Total number of pixels evaluated
        }
        row_acc.next_row(); // Advance to next row
      }
      // Done accumulating valid pixels, now compute the means
      double errorSq = m_max_mean_diffSq + 1.0; // Default value ensures error if there are no valid pixels
      if (matched > 0){
        meanX = meanX / static_cast<double>(matched);
        meanY = meanY / static_cast<double>(matched);
        
        // Compute squared difference of this pixel from the mean disparity
        double thisX   = (*acc)[0];
        double thisY   = (*acc)[1];
        errorSq = (thisX - meanX)*(thisX - meanX) + (thisY - meanY)*(thisY - meanY);
      }
      
      if (errorSq > m_max_mean_diffSq){ // Reject pixels too far from the mean
        m_state->rejected_points++;
        return typename PixelAccessorT::pixel_type();  // Return invalid pixel
      }
    
      return *acc; // Return valid pixel
    }
  }; // End class RmOutliersUsingMeanFunc

  // Useful routine for printing how many points have been rejected
  // using a particular RmOutliersFunc.
  template <class PixelT>
  inline std::ostream&
  operator<<(std::ostream& os, RmOutliersUsingMeanFunc<PixelT> const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices ("
       << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingMeanFunc<typename ViewT::pixel_type> >
  rm_outliers_using_mean(ImageViewBase<ViewT> const& disparity_map,
                         int32 half_h_kernel, int32 half_v_kernel,
                         double max_mean_diff) {
    typedef RmOutliersUsingMeanFunc<typename ViewT::pixel_type> func_type;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, func_type>
      view_type;
    return view_type(edge_extend(disparity_map.impl(), ConstantEdgeExtension()),
                     func_type(half_h_kernel, half_v_kernel,
                               max_mean_diff));
  }

  template <class ViewT>
  inline UnaryPerPixelAccessorView< UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingMeanFunc<typename ViewT::pixel_type> >, RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >
  disparity_cleanup_using_mean(ImageViewBase<ViewT> const& disparity_map,
                               int32 h_half_kernel, int32 v_half_kernel,
                               double max_mean_diff){
    // Rm outliers first using user specified parameters, and then
    // using a heuristic that isolates single pixel outliers.
    typedef RmOutliersUsingThreshFunc<typename ViewT::pixel_type> func_type_thresh;
    typedef RmOutliersUsingMeanFunc<typename ViewT::pixel_type> func_type_mean;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>,
      func_type_mean > inner_type;
    typedef UnaryPerPixelAccessorView<inner_type, func_type_thresh> outer_type;
  
    return outer_type(rm_outliers_using_mean(disparity_map.impl(),
                                             h_half_kernel, v_half_kernel,
                                             max_mean_diff),
                      // Constants set to find only very isolated pixels
                      func_type_thresh( 1, 1, 3.0, 0.2 ) );
  }

  // Method 3: Use the standard deviation of the kernel pixels to
  // determine valid range

  //  rm_outliers_using_stddev()
  //
  // Replacement for old erosion based filter.
  template <class PixelT>
  class RmOutliersUsingStdDev : public ReturnFixedType<PixelT>{

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelAccessorView.
    // By maintaining a smart pointer to this small status class, we
    // can change state that is shared between any copies of the
    // RmOutliersFunc object and the original.
    struct RmOutliersState {
      int32 rejected_points, total_points;
    };

    int32 m_half_h_kernel, m_half_v_kernel;
    double m_pixel_threshold;
    double m_rejection_threshold;
    boost::shared_ptr<RmOutliersState> m_state;

  public:
  
    RmOutliersUsingStdDev(int32 half_h_kernel, int32 half_v_kernel, double pixel_threshold, double rejection_threshold) :
      m_half_h_kernel(half_h_kernel), m_half_v_kernel(half_v_kernel),
      m_pixel_threshold(pixel_threshold), m_rejection_threshold(rejection_threshold),
      m_state( new RmOutliersState() ) {
      m_state->rejected_points = m_state->total_points = 0;

      VW_ASSERT(half_h_kernel > 0 && half_v_kernel > 0,
                ArgumentErr() << "RmOutliersFunc: half kernel sizes must be non-zero.");
    }

    int32 half_h_kernel() const { return m_half_h_kernel; }
    int32 half_v_kernel() const { return m_half_v_kernel; }
    double rejection_threshold() const { return m_rejection_threshold; }
    double pixel_threshold() const { return m_pixel_threshold; }
    int32 rejected_points() const { return m_state->rejected_points; }
    int32 total_points() const { return m_state->total_points; }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_half_h_kernel, -m_half_v_kernel),
                                             Vector2i(m_half_h_kernel, m_half_v_kernel)); }

    // This function will be called for each pixel in the image to be filtered
    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const 
    {
      m_state->total_points++; // Accumulate number of points operated on

      // Quit immediately if the current point is already invalid
      if (!is_valid(*acc))
        return *acc;

      // Allocate storage for recording pixel values
      const size_t numPixels = (2*m_half_v_kernel + 1) * (2*m_half_v_kernel + 1);
      std::vector<double> xVals(numPixels), yVals(numPixels);

      size_t matched = 0, total = 0;
      double meanX=0, meanY=0;
      PixelAccessorT row_acc = acc;
      row_acc.advance(-m_half_h_kernel,-m_half_v_kernel); // Move to top left kernel
      for(int32 yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk) 
        {
          PixelAccessorT col_acc = row_acc; // Start column iterator at left of current row
          for(int32 xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk) 
            {
              if(is_valid(*col_acc))
                {
                  // Record the X and Y values
                  xVals[matched] = (*col_acc)[0];
                  yVals[matched] = (*col_acc)[1];           
                  meanX += (*col_acc)[0];
                  meanY += (*col_acc)[1];
                  matched++; // Total number of valid pixels
                }
              col_acc.next_col(); // Advance to next column
              total++; // Total number of pixels evaluated
            }
          row_acc.next_row(); // Advance to next row
        }
      if (matched == 0) // Reject pixel if there are no valid neighbors
        {
          m_state->rejected_points++;
          return typename PixelAccessorT::pixel_type();  // Return invalid pixel
        }
      // Done accumulating valid pixels, compute means
      meanX = meanX / static_cast<double>(matched);
      meanY = meanY / static_cast<double>(matched);
      
      // Now go through the pixels again to compute the standard deviation
      // - Since we recorded the pixel values in vectors this is easier
      double sumDiffX=0, sumDiffY=0;
      for (size_t i=0; i<matched; ++i)
        {
          double diffX  = xVals[i] - meanX;
          double diffY  = yVals[i] - meanY;
          //double diffSq = diffX*diffX + diffY*diffY;
          sumDiffX += diffX*diffX;
          sumDiffY += diffY*diffY;
        }
      // Compute standard deviation but enforce minimum value
      double stdDevX = sqrt(sumDiffX / static_cast<double>(matched));
      double stdDevY = sqrt(sumDiffY / static_cast<double>(matched));
      if (stdDevX < m_rejection_threshold)
        stdDevX = m_rejection_threshold;
      if (stdDevY < m_rejection_threshold)
        stdDevY = m_rejection_threshold;
      
      // Compute difference of this pixel from the mean disparity
      double thisX = (*acc)[0];
      double thisY = (*acc)[1];
      double errorX = std::abs(thisX - meanX);
      double errorY = std::abs(thisY - meanY);
            
      if ((errorX > m_pixel_threshold*stdDevX) || 
          (errorY > m_pixel_threshold*stdDevY)   ){
        m_state->rejected_points++;
        return typename PixelAccessorT::pixel_type();  // Return invalid pixel
      }
    
      return *acc; // Return valid pixel
    }
  }; // End class RmOutliersUsingStdDev

  // Useful routine for printing how many points have been rejected
  // using a particular RmOutliersFunc.
  template <class PixelT>
  inline std::ostream&
  operator<<(std::ostream& os, RmOutliersUsingStdDev<PixelT> const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices ("
       << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingStdDev<typename ViewT::pixel_type> >
  rm_outliers_using_stddev(ImageViewBase<ViewT> const& disparity_map,
                           int32 half_h_kernel, int32 half_v_kernel,
                           double pixel_threshold,
                           double rejection_threshold) {
    typedef RmOutliersUsingStdDev<typename ViewT::pixel_type> func_type;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, func_type > view_type;
    return view_type(edge_extend(disparity_map.impl(), ConstantEdgeExtension()),
                     func_type(half_h_kernel, half_v_kernel,
                               pixel_threshold, rejection_threshold));
  }

  template <class ViewT>
  inline UnaryPerPixelAccessorView< UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingStdDev<typename ViewT::pixel_type> >, RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >
  disparity_cleanup_using_stddev(ImageViewBase<ViewT> const& disparity_map,
                                 int32 h_half_kernel, int32 v_half_kernel,
                                 double pixel_threshold, double rejection_threshold){
    // Rm outliers first using user specified parameters, and then
    // using a heuristic that isolates single pixel outliers.
    typedef RmOutliersUsingThreshFunc<typename ViewT::pixel_type> func_type_thresh;
    typedef RmOutliersUsingStdDev<typename ViewT::pixel_type> func_type_stddev;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>,
      func_type_stddev > inner_type;
    typedef UnaryPerPixelAccessorView<inner_type,
      func_type_thresh > outer_type;
    return outer_type(rm_outliers_using_stddev(disparity_map.impl(),
                                               h_half_kernel, v_half_kernel,
                                               pixel_threshold, rejection_threshold),
                      func_type_thresh( 1, 1, 3.0, 0.2 ) ); // Constants set to find only very isolated pixels
  }

  // Method 4: Fit a plane to the other pixels and see how well the
  // test pixel fits the plane.
    
  bool fitPlaneToPoints(const std::vector<Vector3> &points, Vector3 &planeDesc);
  double pointToPlaneDist(const Vector3 &point, const Vector3 &planeDesc);
  bool checkPointToPlaneFit(const std::vector<Vector3> &points,
                            const Vector3 &planeDesc,
                            double &meanError, double &stdDevError);
    
  //  rm_outliers_using_plane()
  template <class PixelT>
  class RmOutliersUsingPlane : public ReturnFixedType<PixelT>{

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelAccessorView.
    // By maintaining a smart pointer to this small status class, we
    // can change state that is shared between any copies of the
    // RmOutliersFunc object and the original.
    struct RmOutliersState {
      int32 rejected_points, total_points;
    };

    int32 m_half_h_kernel, m_half_v_kernel;
    double m_pixel_threshold;
    double m_rejection_threshold;
    boost::shared_ptr<RmOutliersState> m_state;

  public:
  
    RmOutliersUsingPlane(int32 half_h_kernel, int32 half_v_kernel, double pixel_threshold, double rejection_threshold) :
      m_half_h_kernel(half_h_kernel), m_half_v_kernel(half_v_kernel),
      m_pixel_threshold(pixel_threshold), m_rejection_threshold(rejection_threshold),
      m_state( new RmOutliersState() ) {
      m_state->rejected_points = m_state->total_points = 0;

      VW_ASSERT(half_h_kernel > 0 && half_v_kernel > 0,
                ArgumentErr() << "RmOutliersFunc: half kernel sizes must be non-zero.");
    }

    int32 half_h_kernel() const { return m_half_h_kernel; }
    int32 half_v_kernel() const { return m_half_v_kernel; }
    double rejection_threshold() const { return m_rejection_threshold; }
    double pixel_threshold() const { return m_pixel_threshold; }
    int32 rejected_points() const { return m_state->rejected_points; }
    int32 total_points() const { return m_state->total_points; }

    BBox2i work_area() const { return BBox2i(Vector2i(-m_half_h_kernel, -m_half_v_kernel),
                                             Vector2i(m_half_h_kernel, m_half_v_kernel)); }

    // This function will be called for each pixel in the image to be filtered
    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const 
    {
      m_state->total_points++; // Accumulate number of points operated on

      // Quit immediately if the current point is already invalid
      if (!is_valid(*acc))
        return *acc;

      // Allocate storage for recording pixel values
      const size_t numPixels = (2*m_half_v_kernel + 1) * (2*m_half_v_kernel + 1);
      std::vector<double> xVals(numPixels), yVals(numPixels);

      // Record all valid points as x/y/z pairs (one set for dX, one set for dY)
      size_t matched = 0, total = 0;
      std::vector<Vector3> xVec, yVec;
      xVec.reserve(numPixels);
      yVec.reserve(numPixels);
      PixelAccessorT row_acc = acc;
      row_acc.advance(-m_half_h_kernel,-m_half_v_kernel); // Move to top left kernel
      for(int32 yk = -m_half_v_kernel; yk <= m_half_v_kernel; ++yk) 
        {
          PixelAccessorT col_acc = row_acc; // Start column iterator at left of current row
          for(int32 xk = -m_half_h_kernel; xk <= m_half_h_kernel; ++xk) 
            {
              if(is_valid(*col_acc))
                {
                  // Record the X and Y values
                  xVec.push_back(Vector3(xk, yk, (*col_acc)[0]));
                  yVec.push_back(Vector3(xk, yk, (*col_acc)[1]));
                  matched++; // Total number of valid pixels
                }
              col_acc.next_col(); // Advance to next column
              total++; // Total number of pixels evaluated
            }
          row_acc.next_row(); // Advance to next row
        }
      if (matched == 0) // Reject pixel if there are no valid neighbors
        {
          m_state->rejected_points++;
          return typename PixelAccessorT::pixel_type();  // Return invalid pixel
        }

      //TODO: eliminate outliers here?

      // Compute a best fit plane for the dX and dY values
      Vector3 xPlane, yPlane;
      try
        {
          fitPlaneToPoints(xVec, xPlane);
          fitPlaneToPoints(yVec, yPlane);
        }
      catch(...) // Failed to solve, probably because points were in a line
        {
          //TODO: Have a fallback check (maybe a line fit) in this case!
          return *acc; // Return valid pixel        
        }

      // Now determine the quality of the plane fits
      double xMean, yMean, xStdDev, yStdDev;
      checkPointToPlaneFit(xVec, xPlane, xMean, xStdDev);
      checkPointToPlaneFit(yVec, yPlane, yMean, yStdDev);
      
      // Enforce minimum standard deviation value
      if (xStdDev < m_rejection_threshold)
        xStdDev = m_rejection_threshold;
      if (yStdDev < m_rejection_threshold)
        yStdDev = m_rejection_threshold;
      
      // Determine if the test pixel is too far from either plane
      Vector3 thisLocX(0,0,(*acc)[0]);
      Vector3 thisLocY(0,0,(*acc)[1]);
      double errorX = pointToPlaneDist(thisLocX, xPlane);
      double errorY = pointToPlaneDist(thisLocY, yPlane);
            
      // Reject pixels greater than N deviations from the mean
      if ( (errorX > m_pixel_threshold*xStdDev) ||
           (errorY > m_pixel_threshold*yStdDev) ){
        m_state->rejected_points++;
        return typename PixelAccessorT::pixel_type();  // Return invalid pixel
      }
    
      return *acc; // Return valid pixel
    }
  }; // End class RmOutliersUsingStdDev

  // Useful routine for printing how many points have been rejected
  // using a particular RmOutliersFunc.
  template <class PixelT>
  inline std::ostream&
  operator<<(std::ostream& os, RmOutliersUsingPlane<PixelT> const& u) {
    os << "\tKernel: [ " << u.half_h_kernel()*2 << ", " << u.half_v_kernel()*2 << "]\n";
    os << "   Rejected " << u.rejected_points() << "/" << u.total_points() << " vertices ("
       << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingPlane<typename ViewT::pixel_type> >
  rm_outliers_using_plane(ImageViewBase<ViewT> const& disparity_map,
                          int32 half_h_kernel, int32 half_v_kernel,
                          double pixel_threshold,
                          double rejection_threshold) {
    typedef RmOutliersUsingPlane<typename ViewT::pixel_type> func_type;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, func_type > view_type;
    return view_type(edge_extend(disparity_map.impl(), ConstantEdgeExtension()),
                     func_type(half_h_kernel, half_v_kernel,
                               pixel_threshold, rejection_threshold));
  }

  template <class ViewT>
  inline UnaryPerPixelAccessorView< UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, RmOutliersUsingPlane<typename ViewT::pixel_type> >, RmOutliersUsingThreshFunc<typename ViewT::pixel_type> >
  disparity_clean_using_plane(ImageViewBase<ViewT> const& disparity_map,
                              int32 h_half_kernel, int32 v_half_kernel,
                              double pixel_threshold, double rejection_threshold) {
    // Rm outliers first using user specified parameters, and then
    // using a heuristic that isolates single pixel outliers.
    typedef RmOutliersUsingThreshFunc<typename ViewT::pixel_type> func_type_thresh;
    typedef RmOutliersUsingPlane<typename ViewT::pixel_type> func_type_plane;
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>,
      func_type_plane > inner_type;
    typedef UnaryPerPixelAccessorView<inner_type,
      func_type_thresh > outer_type;
    return outer_type(rm_outliers_using_plane(disparity_map.impl(),
                                              h_half_kernel, v_half_kernel,
                                              pixel_threshold, rejection_threshold),
                      // Constants set to find only very isolated pixels
                      func_type_thresh( 1, 1, 3.0, 0.2 ) );
  }

  //  std_dev_image()
  //
  /// Remove pixels from the disparity map that correspond to low
  /// contrast pixels in the original image.
  class StdDevImageFunc : public UnaryReturnTemplateType<PixelTypeFromPixelAccessor>
  {
    int32 m_kernel_width, m_kernel_height;

  public:
    StdDevImageFunc(int32 kernel_width, int32 kernel_height);

    BBox2i work_area() const;

    template <class PixelAccessorT>
    typename PixelAccessorT::pixel_type operator() (PixelAccessorT const& acc) const {
      typedef typename PixelAccessorT::pixel_type pixel_type;

      // First pass, compute the mean.
      pixel_type sum = 0;
      PixelAccessorT row_acc = acc;
      row_acc.advance(-m_kernel_width/2,-m_kernel_height/2);
      for(int32 yk = -m_kernel_height/2; yk <= m_kernel_height/2; ++yk) {
        PixelAccessorT col_acc = row_acc;
        for(int32 xk = -m_kernel_width/2; xk <= m_kernel_width/2; ++xk) {
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
      for(int32 yk = -m_kernel_height/2; yk <= m_kernel_height/2; ++yk) {
        PixelAccessorT col_acc = row_acc;
        for(int32 xk = -m_kernel_width/2; xk <= m_kernel_width/2; ++xk) {
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
                int32 kernel_width, int32 kernel_height,
                EdgeT edge) {
    typedef UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, StdDevImageFunc> view_type;
    return view_type(edge_extend(image.impl(), edge),
                     StdDevImageFunc (kernel_width, kernel_height));
  }
  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, StdDevImageFunc>
  std_dev_image(ImageViewBase<ViewT> const& image,
                int32 kernel_width, int32 kernel_height) {
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
  BinaryPerPixelView<ViewT, PerPixelIndexView<VectorIndexFunctor>, TransformDisparitiesFunc<TransformT, typename ViewT::pixel_type> >
  transform_disparities(ImageViewBase<ViewT> const& disparity_map, TransformT const& transform) {

    // Note: We use the PixelIndexView and Binary per pixel filter
    // idiom her to pass the location (in pixel coordinates) into the
    // functor along with the pixel value at that location.
    typedef TransformDisparitiesFunc<TransformT, typename ViewT::pixel_type> func_type;
    typedef BinaryPerPixelView<ViewT, PerPixelIndexView<VectorIndexFunctor>, func_type > view_type;
    return view_type(disparity_map.impl(),pixel_index_view(disparity_map),
                     func_type(transform));
  }



  // Method 5: Use statistics based on a specified region.
  // - This method has some interactions with the tiling system that can result in
  //   large-ish valid portions of the disparity map being flagged as invalid.
  // - More tweaking is required before this can be generally used.

  /// Simple functor to accumulate two CDF functions of disparity data.
  /// - WARNING: The CDF class only works correctly on images of a certain size.
  template <class pixel_type>
  struct DisparityCdfFunctor {
  private:
    size_t count;
    vw::math::CDFAccumulator<float> xCdf, yCdf;
  public:
    DisparityCdfFunctor() : count(0) {}
    
    inline void operator()(pixel_type const& pix) {
      if (!is_valid(pix))
        return; // Skip invalid input pixels
      xCdf(static_cast<float>(pix[0]));
      yCdf(static_cast<float>(pix[1]));
      ++count;
    }
    // The get functions protect against seg faults when () was not called.
    // - This protection belongs in the CDFAccumulator class!
    float getQuantileX(float quantile) const {
      vw_out() << "Quantile count: " << count << std::endl;
      if (count == 0)
        return 0;
      return xCdf.quantile(quantile);
    }
    float getQuantileY(float quantile) const {
      if (count == 0)
        return 0;
      return yCdf.quantile(quantile);
    }
  };

  /// Simple functor to compare disparity values to a pair of thresholds.
  /// - Invalidates the pixel if either value is over the given threshold.
  /// - Could generalize this if desired.
  template <typename T>
  struct DisparityThresholdFunctor {
  private:
    float xCutoff, yCutoff;
  public:
    typedef T result_type;
    
    DisparityThresholdFunctor(float x_cutoff, float y_cutoff) : xCutoff(x_cutoff), yCutoff(y_cutoff) {}
    
    /// Return a copy of the input pixel but invalidate it if it fails either threshold test.
    inline result_type operator()(T const& pix) const {
      T out = pix;
      if ((static_cast<float>(pix[0]) > xCutoff) || (static_cast<float>(pix[1]) > yCutoff))
        invalidate(out);
      return out;
    }

  };

  /// Wrapper to apply the DisparityCdfFunctor plus threshold functor to an entire image.
  template <class ViewT>
  UnaryPerPixelView<ViewT, DisparityThresholdFunctor<typename ViewT::pixel_type> >
  rm_outliers_using_quantiles(ImageViewBase<ViewT> const& disparity_map,
                              double quantile, double multiple) {
                           
    // Compute statistics on the entire input image using DisparityCdfFunctor
    DisparityCdfFunctor<typename ViewT::pixel_type> cdf_functor;
    for_each_pixel(disparity_map, cdf_functor); // Apply the functor to each pixel in the image
    double max_val_x = multiple * cdf_functor.getQuantileX(quantile);
    double max_val_y = multiple * cdf_functor.getQuantileY(quantile);
    //vw_out() << "Computed max value X: " << max_val_x << std::endl;
    //vw_out() << "Computed max value Y: " << max_val_y << std::endl;

    // Use UnaryPerPixelView to apply the threshold functor to every image pixel 
    typedef DisparityThresholdFunctor<typename ViewT::pixel_type> disp_thresh_functor;
    typedef UnaryPerPixelView<ViewT, disp_thresh_functor > view_type;
    return view_type(disparity_map.impl(), disp_thresh_functor(max_val_x, max_val_y));
  }
  
template <class ViewT>
  inline UnaryPerPixelView<ImageView<typename ViewT::pixel_type>, DisparityThresholdFunctor<typename ViewT::pixel_type> > 
  triple_disparity_cleanup(ImageViewBase<ViewT> const& disparity_map,
                           int32  h_half_kernel,   int32  v_half_kernel,
                           double pixel_threshold, double rejection_threshold,
                           double quantile,        double multiple) {

    // Perform speckle filter first and rasterize results so don't need to make an extra pass through them.
    ImageView<typename ViewT::pixel_type> temp_view = 
          disparity_cleanup_using_thresh(disparity_map, h_half_kernel, v_half_kernel,
                                         pixel_threshold, rejection_threshold);

    // Now perform the quantile based filtering.                                         
    return rm_outliers_using_quantiles(temp_view, quantile, multiple); 
  }

template <class ViewT>
  inline UnaryPerPixelView<ImageView<typename ViewT::pixel_type>, DisparityThresholdFunctor<typename ViewT::pixel_type> > 
  double_disparity_cleanup(ImageViewBase<ViewT> const& disparity_map,
                           int32  h_half_kernel,   int32  v_half_kernel,
                           double pixel_threshold, double rejection_threshold,
                           double quantile,        double multiple) {

    // Perform speckle filter first and rasterize results so don't need to make an extra pass through them.
    ImageView<typename ViewT::pixel_type> temp_view = 
          rm_outliers_using_thresh(disparity_map, h_half_kernel, v_half_kernel,
                                         pixel_threshold, rejection_threshold);

    // Now perform the quantile based filtering.                                         
    return rm_outliers_using_quantiles(temp_view, quantile, multiple); 
  }




  // DisparityTransform image transform functor
  //
  // Used to transform an image by using a disparity map
  // Given Left I. + Disparity I. = Right Image :
  // transform( right_image, DisparityTransform(DisparityI)) will
  // project the right image into the perspective of the left.
  class DisparityTransform : public TransformBase<DisparityTransform> {
    typedef ImageViewRef<PixelMask<Vector2f> > inner_view;
    typedef ZeroEdgeExtension EdgeT;
    typedef EdgeExtensionView<inner_view,EdgeT> edge_extend_view;
    typedef NearestPixelInterpolation InterpT;
    InterpolationView<edge_extend_view, InterpT> m_offset_image;
  public:
    template <class DisparityT>
    DisparityTransform(ImageViewBase<DisparityT> const& offset_image) :
      m_offset_image( interpolate(inner_view(offset_image.impl()),InterpT(),EdgeT()) ) {}

    inline Vector2 reverse(const Vector2 &p ) const {
      PixelMask<Vector2f> offset = m_offset_image(p.x(),p.y());
      if ( !is_valid(offset) )
        return Vector2(-1,p.y());
      return p + offset.child();
    }
  };


  // Given a homography transform restricted to a subregion, adjust
  // accordingly the disparity between left and right images.
  template <class DispT>
  ImageView<typename DispT::pixel_type>
  transform_disparities(bool do_round, BBox2i subregion,
                        Matrix<double> const& T, DispT const& disparity){

    VW_ASSERT(subregion.width() == disparity.cols() &&
              subregion.height() == disparity.rows(),
              ArgumentErr() << "transform_disparities: "
              << "The sizes of subregion and disparity don't match.\n");

    ImageView<typename DispT::pixel_type>
      out_disp(disparity.cols(), disparity.rows());
    HomographyTransform H(T);

    for (int x = 0; x < disparity.cols(); x++){
      for (int y = 0; y < disparity.rows(); y++){

        typename DispT::pixel_type disp = disparity(x, y);
        if (!is_valid(disp)) continue; // output stays invalid

        Vector2 beg  = subregion.min() + Vector2(x, y);
        Vector2 end  = beg + disp.child();
        Vector2 tend = H.forward(end);
        Vector2 diff = tend - beg;
        if (do_round) diff = round(diff);
        out_disp(x, y).child() = diff;
        out_disp(x, y).validate();

      }
    }

    return out_disp;
  }

  /// intersect_mask_and_data(view, mask)
  ///
  /// Intersects 'mask' w/ view. View's data is returned first or mask
  template <class PixelT>
  class IntersectPixelMaskData : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    return_type operator()( PixelT const& value, MaskedPixelT const& mask ) const {
      if ( is_valid(value) )
        return value;
      if ( is_valid(mask) )
        return mask;
      return value;
    }
  };

  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMaskData<typename ViewT::pixel_type> >
  intersect_mask_and_data( ImageViewBase<ViewT> const& view,
                           ImageViewBase<MaskViewT> const& mask_view ) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMaskData<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), mask_view.impl(), IntersectPixelMaskData<typename ViewT::pixel_type>() );
  }

  // Disparity Downsample
  template <class ImageT>
  class DisparitySubsampleView : public ImageViewBase<DisparitySubsampleView<ImageT> > {
    ImageT m_child;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename boost::remove_reference<typename ImageT::result_type>::type result_type;
    typedef ProceduralPixelAccessor<DisparitySubsampleView<ImageT> > pixel_accessor;

    DisparitySubsampleView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return 1 + (m_child.cols()-1)/2; }
    inline int32 rows() const { return 1 + (m_child.rows()-1)/2; }
    inline int32 planes() const { return m_child.planes(); }
    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
      int32 ci = i << 1; int32 cj = j << 1;
      typedef typename AccumulatorType<typename PixelChannelType<result_type>::type>::type cnt_type;
      typedef typename CompoundChannelCast<result_type, cnt_type>::type bff_type;
      bff_type buffer;
      cnt_type count = 0;
      if ( is_valid( m_child(ci,cj,p) ) ) {
        count+=10; buffer += 10*bff_type(m_child(ci,cj,p));
      }
      if ( is_valid( m_child(ci+1,cj,p) ) ) {
        count+=5; buffer += 5*bff_type(m_child(ci+1,cj,p));
      }
      if ( is_valid( m_child(ci,cj+1,p) ) ) {
        count+=5; buffer += 5*bff_type(m_child(ci,cj+1,p));
      }
      if ( is_valid( m_child(ci-1,cj,p) ) ) {
        count+=5; buffer += 5*m_child(ci-1,cj,p);
      }
      if ( is_valid( m_child(ci,cj-1,p) ) ) {
        count+=5; buffer += 5*m_child(ci,cj-1,p);
      }
      if ( is_valid( m_child(ci+1,cj+1,p) ) ) {
        count+=2; buffer += 2*m_child(ci+1,cj+1,p);
      }
      if ( is_valid( m_child(ci-1,cj-1,p) ) ) {
        count+=2; buffer += 2*m_child(ci-1,cj-1,p);
      }
      if ( is_valid( m_child(ci-1,cj+1,p) ) ) {
        count+=2; buffer += 2*m_child(ci-1,cj+1,p);
      }
      if ( is_valid( m_child(ci+1,cj-1,p) ) ) {
        count+=2; buffer += 2*m_child(ci+1,cj-1,p);
      }
      if ( count > 0 ) {
        validate(buffer);
        return buffer / (count*2 );
      }
      return result_type();
    }

    typedef DisparitySubsampleView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i expanded_bbox = 2*bbox;
      expanded_bbox.expand(1);
      return prerasterize_type(m_child.prerasterize(expanded_bbox)); }
    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <class ViewT>
  DisparitySubsampleView<EdgeExtensionView<ViewT,ConstantEdgeExtension> >
  disparity_subsample( ImageViewBase<ViewT> const& view ) {
    return DisparitySubsampleView<EdgeExtensionView<ViewT,ConstantEdgeExtension> >( edge_extend(view.impl()) );
  }

  // Disparity Upsample
  template <class ImageT>
  class DisparityUpsampleView : public ImageViewBase<DisparityUpsampleView<ImageT> > {
    ImageT m_child;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename boost::remove_reference<typename ImageT::result_type>::type result_type;
    typedef ProceduralPixelAccessor<DisparityUpsampleView<ImageT> > pixel_accessor;

    DisparityUpsampleView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.cols()*2; }
    inline int32 rows() const { return m_child.rows()*2; }
    inline int32 planes() const { return m_child.planes(); }
    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
      int32 ci = i >> 1; int32 cj = j >> 1;
      return m_child(ci,cj,p)*2;
    }

    typedef DisparityUpsampleView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type(m_child.prerasterize(bbox)); }
    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <class ViewT>
  DisparityUpsampleView<ViewT>
  disparity_upsample( ImageViewBase<ViewT> const& view ) {
    return DisparityUpsampleView<ViewT>( view.impl() );
  }

}}    // namespace vw::stereo


#endif //  _VWDISPARITYMAP_H_
