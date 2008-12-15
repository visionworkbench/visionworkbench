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

// SURF.h
//
// Classes to perform SURF interest point detection. A lot of SURF
// doesn't fit into the cookie cutter model of IP.
//
// Author:
// Zachary Moratto
// Kansas State University
#ifndef __VW_INTERESTPOINT_SURF_H__
#define __VW_INTERESTPOINT_SURF_H__

#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/Math.h>
#include <vw/Image.h>

#include <vector>
#include <list>

namespace vw {
namespace ip {

  /// Red Herring
  class SURFInterestOperator {
    
    static const float DEFAULT_INTEREST_THRESHOLD = 0.03;
    float m_threshold;
  public:
    SURFInterestOperator() : m_threshold(DEFAULT_INTEREST_THRESHOLD) {}
    SURFInterestOperator( float threshold ) : m_threshold(threshold) {}
   
    /// DO NOT USE, this is only for compatibility with IP code
    template <class ViewT>
    inline ImageViewRef< typename ViewT::pixel_type>
    operator() (ImageViewBase<ViewT> const& source, float scale= 1.0) const {
      vw_out(DebugMessage, "interest_point") << "DO NOT USE SURFInterestOperator it does nothing.\n";
    }

    template <class DataT>
    inline bool threshold (InterestPoint const& pt, DataT const& data) const {
      return false;
    }

    template <class DataT>
    inline bool threshold (DataT const& data) const {
      return data > m_threshold;
    }

    float threshold ( void ) const {
      return m_threshold;
    }
  };

  /// SURFParams
  /// Just a fancier way to share data with the functions
  class SURFParams {
  public:
    SURFParams( float thres = .1, unsigned oct = 4, unsigned scal = 4 ) {
      cols = rows = 0;    // Remember to set this later
      threshold = thres;
      octaves = oct;
      scales = scal;
      lobe_size = 3;
      
      // Building indexes
      indexes = new unsigned[octaves * scales];
      unsigned count = 0;
      for (unsigned o = 0; o < octaves; o++ ) {
	for (unsigned s = 0; s < scales; s++ ) {
	  if ( o == 0 ) {
	    indexes[o*scales + s] = count;
	    count++;
	  } else if ( s == 0 ) {
	    indexes[o*scales] = indexes[o*scales - 3];
	  } else if ( s == 1 ) {
	    indexes[o*scales+1] = indexes[o*scales-1];
	  } else {
	    indexes[o*scales + s] = count;
	    count++;
	  }
	}
      }

    }
    template <class ViewT>
    void setForImage( ImageViewBase<ViewT> const& source ) {
      ImageView<float> image = pixel_cast<PixelGray<float> >(source);

      cols = image.cols();
      rows = image.rows();
    }

    template <class T>
    float calcScale( T const& filter_size ) const {
      return 1.2*float(filter_size)/9.0;
    }

    unsigned convToIndex( unsigned const& octave, unsigned const& scale ) const {
      unsigned idx = octave*scales + scale;
      VW_DEBUG_ASSERT( idx < 16, 
		       vw::ArgumentErr() << "convToIndex: Incorrect octave and scale given" );
      return indexes[idx];
    }
    int cols;
    int rows;
    float threshold;
    unsigned octaves;
    unsigned scales;
    unsigned lobe_size;
  private:
    // This is used for converting octave and scale to an index of the
    // interest data. It is fixed for 4 octaves and 4 scales right now
    unsigned *indexes;
  };

  /// SURFScaleData
  /// This is a wrapper for the interest data to make it easier to
  /// compare different scales to each other
  class SURFScaleData {
  public:
    // Fancy constructor
    SURFScaleData( boost::shared_ptr<vw::ImageView<float> > data,
		   boost::shared_ptr<vw::ImageView<bool> > polarity,
		   unsigned const& filter_size, unsigned const& sampStep) :
    m_data(data), m_polarity(polarity), m_filter_size(filter_size), m_sampling_step(sampStep)
    { 
      // calculating sampling_step_2, (finds x of 2^x = sampling_step)
      unsigned temp = m_sampling_step;
      m_sampling_step_2 = 0;
      while ( (temp & 0x1) == 0 && temp != 0 ) {
	m_sampling_step_2++;
	temp>>=1;
      }

      // Calculating the starting offset
      unsigned s_temp = floor(float(m_filter_size)/2.);
      m_starting_offset = s_temp + (m_sampling_step - (s_temp & ( m_sampling_step - 1)));
    }

    float determinant( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 ){
	return 0.0;
      }
      return m_data->operator()(cov[0],cov[1]);
    }

    bool polarity( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return m_polarity->operator()(cov[0],cov[1]);
    }

    bool isValid( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return true;
    }
    
    // isLessThan
    // - value = value to compare interest image against
    // - x,y   = point to compare value against 
    //           (in the coordinates of the orginal image)
    bool isLessThan( float const& value,
		     int const& x, int const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return value < m_data->operator()(cov[0],cov[1]);
    }

    unsigned filter_size( void ) const {
      return m_filter_size;
    }
  private:
    unsigned m_filter_size;
    unsigned m_starting_offset;
    unsigned m_sampling_step;
    unsigned m_sampling_step_2;
    boost::shared_ptr<vw::ImageView<float> > m_data;
    boost::shared_ptr<vw::ImageView<bool> > m_polarity;

    vw::Vector3i convertCoord( int const& x, int const& y) const {
      if ( (x & 1) || (y & 1) ) // is it odd?
	return vw::Vector3i(0,0,0);
      int tx = x, ty = y;
      tx -= m_starting_offset;
      ty -= m_starting_offset;
      if (tx < 0 || ty < 0) // Have we left our space?
	return vw::Vector3i(0,0,0);
      if (tx & (m_sampling_step - 1) || ty & (m_sampling_step - 1) ) // mod not equal zero
	return vw::Vector3i(0,0,0);
      tx = tx >> m_sampling_step_2;
      ty = ty >> m_sampling_step_2;
      if ( tx >= m_data->cols() || ty >= m_data->rows() )
	return vw::Vector3i(0,0,0);
      return vw::Vector3i(tx,ty,1);
    }
  };

  /// Creates a interest scale for SURF
  SURFScaleData SURFProcessScale( vw::ImageView<double> const& integral, 
				  unsigned const& filter_size,
				  unsigned const& sampling_step,
				  SURFParams const& params );

  /// Finds the maximas in our data
  std::list<InterestPoint> SURFMaximaDetection( std::vector<SURFScaleData> const & scaleData, 
						SURFParams const& params );
  
  /// Writes SURF debug interest images
  void SURFWriteDebugImages ( std::vector<SURFScaleData> const& scaleData,
			      SURFParams const& params );

  /// SURF Subpixel refinement
  /// This gets floating point precision for the interest points, will
  /// also throw at some points if they move too far.
  void SURFSubpixelRefinement( std::list<vw::ip::InterestPoint>& ip,
			       std::vector<SURFScaleData> const& scaleData,
			       SURFParams const& params );

  // SURF Gradient 3D
  // This calculates the gradient, hmm... it negative
  Vector3 SURFGradient3D( std::vector<SURFScaleData> const& scaleData,
			  vw::ip::InterestPoint const& ip,
			  SURFParams const& params );

  // SURF Hessian 3D
  // This calculates the hessian matrix for sub pixel correlation
  Matrix3x3 SURFHessian3D( std::vector<SURFScaleData> const& scaleData, 
			   vw::ip::InterestPoint const& ip, 
			   SURFParams const& params );

  // SURF Orientation
  // This calculates the orientation of a feature
  float SURFOrientation( vw::ImageView<double> const&, int const&, int const&,
			 float const& );  

  /// This class actually performs all the work. The SURF operator is
  /// actually just a red herring that is required to match the IP module's
  /// frame.
  template <class InterestT>
  class SURFInterestPointDetector : public InterestDetectorBase<SURFInterestPointDetector<InterestT> >
  {

    SURFInterestPointDetector(SURFInterestPointDetector<InterestT> const& copy) {}

  public:
    static const int IP_DEFAULT_SCALES = 4; // Don't play with
    static const int IP_DEFAULT_OCTAVES = 4;

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwies, the max_points most "interesting" points are
    /// returned.
    SURFInterestPointDetector(int max_points = 0)
      : m_max_points(max_points) {
      m_params = SURFParams();
    }

    SURFInterestPointDetector(InterestT const& interest, int max_points = 0) 
      : m_max_points(max_points) {
      vw_out(DebugMessage, "interest_point") << "SURFInterestPointDetector rejected outside interest operator.\n";
      m_params = SURFParams( interest.threshold(), IP_DEFAULT_OCTAVES, IP_DEFAULT_SCALES );
    }

    SURFInterestPointDetector(InterestT const& interest, int scales, int octaves, int max_points = 0)
      : m_max_points(max_points) {
      vw_out(DebugMessage, "interest_point") << "SURFInterestPointDetector rejected outside interest operator. Also, it is not recommended to change octave and scale values.\n";
      m_params = SURFParams( interest.threshold(), octaves, scales );
    }

    /// Process Image ////////////////////////////
    /// Detect interest points in the source image.
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image) {

      // Timing
      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      // Set SURF Params;
      m_params.setForImage( image );
      
      /*
      // DEBUG TEST
      ImageView<int> test_image(10, 10);
      for (unsigned x = 0; x < 10; x++ )
	for (unsigned y = 0; y < 10; y++ )
	  test_image(x,y) = 2*(x+3) + 3*(y+1);
      ImageView<double> test_integral = IntegralImage( test_image );
      float Dxx = XSecondDerivative( test_integral, 4, 4, 9 );
      float Dyy = YSecondDerivative( test_integral, 4, 4, 9 );
      */

      // Create Integral Image
      vw_out(InfoMessage,"interest_point") << "\tBuilding Integral Image.\n";
      ImageView<double> integral = IntegralImage( image );
      vw_out(DebugMessage,"interest_point") << "\tdone ( integral size " << integral.cols() << "x" << integral.rows() << " ).\n";

      // Create Interest Scales
      vw_out(InfoMessage,"interest_point") << "\tCalculating Interest Data.\n";
      std::vector<SURFScaleData> interest_scales;
      unsigned filter_size = 3;
      for ( unsigned octave = 0; octave < m_params.octaves; octave++ ) {
	unsigned sampling_step = 0x2 << octave;
	unsigned iterations = 2;
	if ( octave == 0 )
	  iterations = 4;
       
	for ( unsigned i = 0; i < iterations; i++ ) {
	  filter_size += 6 * ( octave == 0 ? 1 : 0x2 << ( octave -1 ) );
	  
	  if ( filter_size + 2*sampling_step > integral.cols() ||
	       filter_size + 2*sampling_step > integral.rows() ) {
	    m_params.octaves--;
	    break;
	  }

	  interest_scales.push_back( SURFProcessScale( integral,
						       filter_size,
						       sampling_step,
						       m_params ) );
	}
      }
      
      // Write Debug Images ?
      // SURFWriteDebugImages( interest_scales, m_params );
      
      // Locate Maximas
      vw_out(InfoMessage,"interest_point") << "\tPerforming non-maximal suppression.\n";
      std::list<InterestPoint> ip = SURFMaximaDetection( interest_scales,
							 m_params );

      vw_out(InfoMessage,"interest_point") << "\t\toutput points: " << ip.size() << std::endl;
      
      // Subpixel refinement
      vw_out(InfoMessage,"interest_point") << "\tPerforming sub-pixel refinement.\n";
      SURFSubpixelRefinement( ip, interest_scales, m_params );
      
      vw_out(InfoMessage,"interest_point") << "\t\toutput points: " << ip.size() << std::endl;
      
      // Cull
      if ( m_max_points != 0 ) {
	vw_out(InfoMessage,"interest_point") << "\tCulling data to a " << m_max_points << " interest points.\n";
	int orginal_num_points = ip.size();
	ip.sort();
	if ( m_max_points < ip.size() )
	  ip.resize( m_max_points );
	vw_out(DebugMessage,"interest_point") << "\tdone (removed " << orginal_num_points - ip.size() << " interest points, " << ip.size() << " remaining.).\n";
      }
      
      // Assign orientations
      vw_out(InfoMessage,"interest_point") << "\tCalculating orientation.\n";
      for (std::list<InterestPoint>::iterator point = ip.begin();
	   point != ip.end(); ++point) {
	(*point).orientation = SURFOrientation( integral,
						(*point).ix, (*point).iy,
						(*point).scale );
      }
      
      delete total;
      return ip;
    }


  protected:
    int m_max_points;
    SURFParams m_params;
  };

  // SURF Descriptor 64 bit
  struct SURFDescriptorGenerator : public DescriptorGeneratorBase<SURFDescriptorGenerator> {
    

    Matrix<float,20,20> gaussian_weight;
    bool m_extended;

    // Constructor
    SURFDescriptorGenerator( bool extended = false ) {
      m_extended = extended;
      // building gaussian weight meu = 9.5 (center of our gaussian)
      for (char x = 0; x < 20; x++ ) {
	for (char y = 0; y < 20; y++ ) {
	  float dist = (float(x) - 9.5f)*(float(x) - 9.5 ) + (float(y) - 9.5)*(float(y) - 9.5 );
	  
	  gaussian_weight(x,y) = exp( -dist/21.78 );
	}
      }
    }

    // Actually Descriptor arithmetic
    template <class ViewT>
    Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) const {
      
      Vector<float> result(64);
      if (m_extended)
	result.set_size(128);
      Matrix<float,20,20> h_response;
      Matrix<float,20,20> v_response;

      // Building responses
      for ( char x = 0; x < 20; x++ ) {
	for ( char y = 0; y < 20; y++ ) {
	  h_response(x,y) = (support.impl()(x+1,y).v() + support.impl()(x+1,y+1).v() - 
			     support.impl()(x,y).v() - support.impl()(x,y+1).v() )*gaussian_weight(x,y);

	  v_response(x,y) = (support.impl()(x,y+1).v() + support.impl()(x+1,y+1).v() - 
			     support.impl()(x,y).v() - support.impl()(x+1,y).v() )*gaussian_weight(x,y);
	}
      }

      // Building the descriptor
      for ( char x = 0; x < 4; x++ ) {
	for ( char y = 0; y < 4; y++ ) {
	  
	  // Summing responses
	  for ( char ix = 0; ix < 5; ix++ ) {
	    for ( char iy = 0; iy < 5; iy++ ) {

	      if (m_extended) {
		//SURF128

		// Working the DXs
		if ( v_response(x*5+ix, y*5+iy) > 0 ) {
		  result(32*x + 8*y + 0) = h_response(x*5+ix, y*5+iy);
		  result(32*x + 8*y + 1) = fabs(h_response(x*5+ix, y*5+iy));
		} else {
		  result(32*x + 8*y + 2) = h_response(x*5+ix, y*5+iy);
		  result(32*x + 8*y + 3) = fabs(h_response(x*5+ix, y*5+iy));
		}

		// Working the DYs
		if ( h_response(x*5+ix, y*5+iy) > 0 ) {
		  result(32*x + 8*y + 4) = v_response(x*5+ix, y*5+iy);
		  result(32*x + 8*y + 5) = fabs(v_response(x*5+ix, y*5+iy));
		} else {
		  result(32*x + 8*y + 6) = v_response(x*5+ix, y*5+iy);
		  result(32*x + 8*y + 7) = fabs(v_response(x*5+ix, y*5+iy));
		}

	      } else {
		//SURF 64 (The Orginal)

		// Dx
		result(16*x + 4*y) += h_response(x*5+ix, y*5+iy);

		// |Dx|
		result(16*x + 4*y + 1) += fabs(h_response(x*5+ix, y*5+iy));

		// Dy
		result(16*x + 4*y + 2) += v_response(x*5+ix, y*5+iy);
	     
		// |Dy|
		result(16*x + 4*y + 3) += fabs(v_response(x*5+ix, y*5+iy));
	      }
	    }
	  }

	}
      }

      return normalize(result);
    }

    // This builds the window used for the descriptor
    template <class ViewT>
    inline ImageView<typename ViewT::pixel_type> get_support( float x, float y, float scale, float ori,
                                                              ImageViewBase<ViewT> const& source, int size=21 ) {

      float scaling = 1.0f / scale;

      // Output is 21x21 pixels
      return transform(source.impl(),
		       compose(TranslateTransform(10, 10), 
			       compose(ResampleTransform(scaling, scaling),
				       RotateTransform(-ori),
				       TranslateTransform(-x, -y))),
		       21, 21); 
    }

  };
  
}} // namespace vw::ip

#endif // __VW_INTERESTPOINT_SURF_H__
