// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file IntegralDetector.h
///
/// Detectors that operate off of the use of integral images
///

#ifndef __VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__
#define __VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__

#include <vw/Image/Interpolation.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralImage.h>
#include <deque>

namespace vw {
namespace ip {

  /// Getting Orientation. This method of finding orientation is the
  /// same as described in SURF papers.
  template <class IntegralT>
  struct AssignOrientation {
    IntegralT integral;
    AssignOrientation( ImageViewBase<IntegralT> const& image ) : integral(image.impl()) {}

    void operator()( InterestPoint& ip ) {
      typedef std::vector<typename IntegralT::pixel_type> Measures;
      Measures h_response(169);
      Measures v_response(169);
      Measures angle(169);
      int m = 0;

      typedef InterpolationView<EdgeExtensionView< IntegralT, ConstantEdgeExtension>, BilinearInterpolation> WrappedType;
      WrappedType wrapped_integral = interpolate(integral.impl());

      for ( int i = -6; i <= 6; i++ ) {
        for ( int j = -6; j <= 6; j++ ) {
          float distance_2 = i*i+j*j;

          // Check if the location is within the radius of 6
          if ( distance_2 > 36 )
            continue;
          float weight = exp(-distance_2/8)/5.0133;

          Vector2 location = Vector2(ip.ix,ip.iy) +
            ip.scale*Vector2(i,j);
          h_response[m] = weight*HHaarWavelet( wrapped_integral,
                                               location[0],
                                             location[1],
                                               ip.scale*4 );
          v_response[m] = weight*VHaarWavelet( wrapped_integral,
                                               location[0],
                                               location[1],
                                               ip.scale*4 );
          angle[m] = atan2( v_response[m],
                            h_response[m] );
          m++;
        }
      }

      // Fitting a slice to find response
      static const float pi_6 = M_PI/6.0;
      static const float two_pi = 2*M_PI;
      float sumx, sumy;
      float mod, greatest_mod=0;
      float greatest_ori=0;
      for ( float a = 0; a < two_pi; a += 0.5 ) {
        sumx = sumy = 0;
        for ( int idx = 0; idx < m; idx++ ) {
          // Is it in my slice
          if ( (angle[idx] > a - pi_6 && angle[idx] < a + pi_6 ) ||
               (angle[idx] + two_pi > a - pi_6 && angle[idx] + two_pi < a + pi_6 ) ||
               (angle[idx] - two_pi > a - pi_6 && angle[idx] - two_pi < a + pi_6 ) ) {
            sumx += h_response[idx];
            sumy += v_response[idx];
          }
        }

        mod = sumx*sumx+sumy*sumy;
        if ( mod > greatest_mod ) {
          greatest_mod = mod;
          greatest_ori = atan2(sumy,sumx);
        }
      }

      ip.orientation = greatest_ori;
    }
  };

  // IntegralInterestPointDetector
  template <class InterestT>
  class IntegralInterestPointDetector : public InterestDetectorBase<IntegralInterestPointDetector<InterestT> >,
                                        private boost::noncopyable {

  public:
    static const int IP_DEFAULT_SCALES = 8;

    /// Setting max_points = 0 will disable interest point culling.
    IntegralInterestPointDetector( int max_points = 1000 )
      : m_interest(InterestT()), m_scales(IP_DEFAULT_SCALES), m_max_points(max_points) {}

    IntegralInterestPointDetector( InterestT const& interest, int max_points = 1000 )
      : m_interest(interest), m_scales(IP_DEFAULT_SCALES), m_max_points(max_points) {}

    IntegralInterestPointDetector( InterestT const& interest, int scales, int max_points )
      : m_interest(interest), m_scales(scales), m_max_points(max_points) {}

    /// Detect Interest Points in the source image.
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image ) const {
      typedef ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> ImageT;
      typedef ImageInterestData<ImageT,InterestT> DataT;

      Timer total("\t\tTotal elapsed time", DebugMessage, "interest_point");

      // Rendering own standard copy of the image as the passed in view is just a cropview
      ImageT original_image = image.impl();

      // Producing Integral Image
      ImageT integral_image;
      {
        vw_out(DebugMessage, "interest_point") << "\tCreating Integral Image ...";
        Timer t("done, elapsed time", DebugMessage, "interest_point");
        integral_image= IntegralImage( original_image );
      }

      // Creating Scales
      std::deque<DataT> interest_data;
      interest_data.push_back( DataT(original_image, integral_image) );
      interest_data.push_back( DataT(original_image, integral_image) );

      // Priming scales
      InterestPointList new_points;
      {
        vw_out(DebugMessage, "interest_point") << "\tScale 0 ... ";
        Timer t("done, elapsed time", DebugMessage, "interest_point");
        m_interest( interest_data[0], 0 );
      }
      {
        vw_out(DebugMessage, "interest_point") << "\tScale 1 ... ";
        Timer t("done, elapsed time", DebugMessage, "interest_point");
        m_interest( interest_data[1], 1 );
      }
      // Finally processing scales
      for ( int scale = 2; scale < m_scales; scale++ ) {

        interest_data.push_back( DataT(original_image, integral_image) );
        {
          vw_out(DebugMessage, "interest_point") << "\tScale " << scale << " ... ";
          Timer t("done, elapsed time", DebugMessage, "interest_point");
          m_interest( interest_data[2], scale );
        }

        InterestPointList scale_points;

        // Detecting interest points in middle
        int32 cols = original_image.cols() - 2;
        int32 rows = original_image.rows() - 2;
        typedef typename DataT::interest_type::pixel_accessor AccessT;

        AccessT l_row = interest_data[0].interest().origin();
        AccessT m_row = interest_data[1].interest().origin();
        AccessT h_row = interest_data[2].interest().origin();
        l_row.advance(1,1); m_row.advance(1,1); h_row.advance(1,1);
        for ( int32 r=0; r < rows; r++ ) {
          AccessT l_col = l_row;
          AccessT m_col = m_row;
          AccessT h_col = h_row;
          for ( int32 c=0; c < cols; c++ ) {
            if ( is_extrema( l_col, m_col, h_col ) )
              scale_points.push_back(InterestPoint(c+2,r+2,
                                                   m_interest.float_scale(scale-1),
                                                   *m_col) );
            l_col.next_col();
            m_col.next_col();
            h_col.next_col();
          }
          l_row.next_row();
          m_row.next_row();
          h_row.next_row();
        }

        // Thresholding
        threshold(scale_points, interest_data[1], scale-1);

        // Appending to the greater set
        new_points.insert(new_points.end(),
                          scale_points.begin(),
                          scale_points.end());

        // Deleting lowest
        interest_data.pop_front();
      }

      if ( m_max_points > 0 ) { // Cull
        vw_out(DebugMessage, "interest_point") << "\tCulling ...";
        Timer t("elapsed time", DebugMessage, "interest_point");
        int original_num_points = new_points.size();
        new_points.sort();
        if (m_max_points < original_num_points)
          new_points.resize( m_max_points );
        vw_out(DebugMessage, "interest_point") << "     (removed " << original_num_points - new_points.size() << " interest points, " << new_points.size() << " remaining.)\n";
      }

      { // Assign orientations
        vw_out(DebugMessage, "interest_point") << "\tAssigning Orientations... ";
        Timer t("elapsed time", DebugMessage, "interest_point");
        std::for_each( new_points.begin(), new_points.end(),
                       AssignOrientation<ImageT >( integral_image ) );
      }

      return new_points;
    }

  protected:

    InterestT m_interest;
    int m_scales, m_max_points;

    template <class AccessT>
    bool inline is_extrema( AccessT const& low,
                            AccessT const& mid,
                            AccessT const& hi ) const {
      AccessT low_o = low;
      AccessT mid_o = mid;
      AccessT hi_o  = hi;

      if ( *mid_o <= *low_o ) return false;
      if ( *mid_o <= *hi_o  ) return false;

      for ( uint8 step = 0; step < 8; step++ ) {
        if ( step == 0 ) {
          low_o.advance(-1,-1); mid_o.advance(-1,-1);hi_o.advance(-1,-1);
        } else if ( step == 1 || step == 2 ) {
          low_o.next_col(); mid_o.next_col(); hi_o.next_col();
        } else if ( step == 3 || step == 4 ) {
          low_o.next_row(); mid_o.next_row(); hi_o.next_row();
        } else if ( step == 5 || step == 6 ) {
          low_o.prev_col(); mid_o.prev_col(); hi_o.prev_col();
        } else {
          low_o.prev_row(); mid_o.prev_row(); hi_o.prev_row();
        }
        if ( *mid <= *low_o ||
             *mid <= *mid_o ||
             *mid <= *hi_o ) return false;
      }

      return true;
    }

    template <class DataT>
    inline void threshold( InterestPointList& points,
                           DataT const& img_data,
                           int const& scale ) const {
      InterestPointList::iterator pos = points.begin();
      while (pos != points.end()) {
        if (!m_interest.threshold(*pos,
                                  img_data, scale) )
          pos = points.erase(pos);
        else
          pos++;
      }
    }
  };


}} // end vw::ip

#endif//__VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__
