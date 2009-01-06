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

// IntegralImage.h
//
// Provides support for algorithms that evaluate with Integral Images
// such as SURF
#ifndef __VW_INTERESTPOINT_INTEGRALIMAGE_H__
#define __VW_INTERESTPOINT_INTEGRALIMAGE_H__

#include <boost/utility/enable_if.hpp>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Math.h>

namespace vw {
namespace ip {

  /// Creates Integral Image
  template <class ViewT>
  inline ImageView<double> IntegralImage( ImageViewBase<ViewT> const& source ) {
    
    // Allocating space for integral image
    vw::ImageView<double> integral( source.impl().cols()+1, source.impl().rows()+1 );
    
    // Performing cumulative sum in the x direction
    for( signed y = 1; y < integral.rows(); ++y ) {
      integral(0,y) = 0;  // This is only need if the integral is not zero in the beginning
      for( signed x = 1; x < integral.cols(); ++x ) {
	integral(x,y) = integral(x-1,y) + pixel_cast<PixelGray<double> >(source.impl()(x-1,y-1)).v();
      }
    }

    // Performing cumulative sum in the y direction
    for( int x = 1; x < integral.cols(); ++x ) {
      integral(x,0) = 0;
      for ( int y = 1; y < integral.rows(); ++y)
	integral(x,y) += integral(x,y-1);
    }

    return integral;
  }

  /// Integral Block Evaluation
  ///
  /// This is for summing an area of pixels described by integral
  inline double IntegralBlock( ImageView<double> const& integral, 
			       Vector2i const& top_left,
			       Vector2i const& bottom_right ) {
    VW_DEBUG_ASSERT(top_left.x() > bottom_right.x() && top_left.y() > bottom_right.y(), 
		    vw::ArgumentErr() << "Incorrect input for IntegralBlock.\n");
    VW_DEBUG_ASSERT(top_left.x() < (unsigned)integral.cols(), 
		    vw::ArgumentErr() << "x0 out of bounds. "<< integral.cols() <<" : "
		    << top_left << bottom_right << "\n");
    VW_DEBUG_ASSERT(bottom_right.x() < (unsigned)integral.cols(), 
		    vw::ArgumentErr() << "x1 out of bounds. "<< integral.cols() <<" : "
		    << top_left << bottom_right << "\n");
    VW_DEBUG_ASSERT(top_left.y() < (unsigned)integral.rows(), 
		    vw::ArgumentErr() << "y0 out of bounds. "<< integral.rows() <<" : "
		    << top_left << bottom_right << "\n");
    VW_DEBUG_ASSERT(bottom_right.y() < (unsigned)integral.rows(), 
		    vw::ArgumentErr() << "y1 out of bounds. "<< integral.rows() <<" : "
		    << top_left << bottom_right << "\n");
    
    double result;
    result = integral( top_left.x(), top_left.y() );
    result += integral( bottom_right.x(), bottom_right.y() );
    result -= integral( top_left.x(), bottom_right.y() );
    result -= integral( bottom_right.x(), top_left.y() );

    //result /= (top_left.x() - bottom_right.x())*(top_left.y() - bottom_right.y());

    return result;
  }

  /// X First Derivative
  /// - x,y         = location to center the calculation on
  /// - filter_size = size of window for calculation
  inline float XFirstDerivative( ImageView<double> const& integral,
				 int const& x, int const& y,
				 unsigned const& filter_size ) {
    float derivative = 0;
    return derivative;
  }

  /// X Second Derivative
  /// - x,y         = location to center the calculation on
  /// - filter_size = size of window for calculation
  inline float XSecondDerivative( ImageView<double> const& integral,
				  int const& x, int const& y,
				  unsigned const& filter_size ) {
    unsigned lobe = filter_size / 3;
    unsigned half_lobe = floor( float(lobe) / 2.0 );
    float derivative;

    // Adding positive left;
    derivative = IntegralBlock( integral,
				Vector2i( x - lobe - half_lobe, y - lobe + 1 ),
				Vector2i( x - half_lobe, y + lobe ) );

    // Adding negative middle;
    derivative -= 2.0*IntegralBlock( integral,
				     Vector2i( x - half_lobe, y - lobe + 1 ),
				     Vector2i( x + half_lobe + 1, y + lobe ) );

    // Adding positive right;
    derivative += IntegralBlock( integral,
				 Vector2i( x + half_lobe + 1, y - lobe + 1 ),
				 Vector2i( x + half_lobe + lobe + 1, y + lobe ) );
    
    derivative /= filter_size*filter_size;

    return derivative;
  }

  /// Y First Derivative
  /// - x,y         = location to center the calculation on
  /// - filter_size = size of window for calculation
  inline float YFirstDerivative( ImageView<double> const& integral,
				 int const& x, int const& y,
				 unsigned const& filter_size ) {
    float derivative = 0;

    return derivative;
  }

  /// Y Second Derivative
  /// - x,y         = location to center the calculation on
  /// - filter_size = size of window for calculation
  inline float YSecondDerivative( ImageView<double> const& integral,
				  int const& x, int const& y,
				  unsigned const& filter_size ) {
    unsigned lobe = filter_size / 3;
    unsigned half_lobe = floor( float(lobe) / 2.0 );
    float derivative;

    // Adding positive top;
    derivative = IntegralBlock( integral,
				Vector2i( x - lobe + 1, y - lobe - half_lobe ),
				Vector2i( x + lobe, y - half_lobe ) );

    // Adding negative middle;
    derivative -= 2.0*IntegralBlock( integral,
				     Vector2i( x - lobe + 1, y - half_lobe ),
				     Vector2i( x + lobe, y + half_lobe + 1 ) );

    // Adding positive bottom;
    derivative += IntegralBlock( integral,
				 Vector2i( x - lobe + 1, y + half_lobe + 1 ),
				 Vector2i( x + lobe, y + half_lobe + lobe + 1 ) );

    derivative /= filter_size*filter_size;

    return derivative;
  }

  /// XY Derivative
  /// - x,y         = location to center the calculation on
  /// - filter_size = size of window for calculation
  inline float XYDerivative( ImageView<double> const& integral,
			     int const& x, int const& y,
			     unsigned const& filter_size ) {

    unsigned lobe = filter_size / 3;
    float derivative;

    // Adding positive top left
    derivative = IntegralBlock( integral,
				Vector2i( x - lobe, y - lobe ),
				Vector2i( x, y ) );
    
    // Adding negative top right
    derivative -= IntegralBlock( integral,
				 Vector2i( x + 1, y - lobe ),
				 Vector2i( x + lobe + 1, y ) );

    // Adding negative bottom left
    derivative -= IntegralBlock( integral,
				 Vector2i( x - lobe, y + 1 ),
				 Vector2i( x, y + lobe + 1 ) );

    // Adding positve bottom right
    derivative += IntegralBlock( integral,
				 Vector2i( x + 1, y + 1 ),
				 Vector2i( x + 1 + lobe, y + 1 + lobe ) );

    derivative /= filter_size*filter_size;

    return derivative;
  }

  // Horizontal Wavelet
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate

  // Note: Filter will be evaluated at a size nearest to a multiple of two
  template <class ViewT, class NumberT>
  typename boost::enable_if<boost::is_integral<NumberT>, float>::type
  inline HHaarWavelet( ImageViewBase<ViewT> const& integral,
			     NumberT const& x, NumberT const& y,
			     float const& size ) {

    float response;
    int half_size = round( size / 2.0);
    int i_size = half_size << 1;
    int top = round( int(y) - size/2);
    int left = round( int(x) - size/2);

    VW_ASSERT(left+i_size < (unsigned)integral.impl().cols(), 
	      vw::ArgumentErr() << "left out of bounds. "<< integral.impl().cols() <<" : "
	      << left+i_size << " [top left] " << top << " " << left << "\n");
    VW_ASSERT(top+i_size < (unsigned)integral.impl().rows(),
	      vw::ArgumentErr() << "top out of bounds. " << integral.impl().rows() <<" : "
	      << top+i_size << "\n");
    VW_ASSERT(left >= 0,
	      vw::ArgumentErr() << "left is to low. " << 0 << " : " << left << "\n");
    VW_ASSERT(top >= 0,
	      vw::ArgumentErr() << "top is to low. " << 0 << " : " << top << "\n");

    response = -integral.impl()(left, top);
    response += 2*integral.impl()(left+half_size, top);
    response -= integral.impl()(left+i_size, top);
    response += integral.impl()(left, top+i_size);
    response -= 2*integral.impl()(left+half_size, top+i_size);
    response += integral.impl()(left+i_size, top+i_size);
    
    return response;
  }

  // Vertical Wavelet
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate
  // Note: Filter will be evaluated at a size nearest to a multiple of two
  template <class ViewT, class NumberT>
  typename boost::enable_if<boost::is_integral<NumberT>, float>::type
  inline VHaarWavelet( ImageViewBase<ViewT> const& integral,
			     NumberT const& x, NumberT const& y,
			     float const& size ) {

    float response;
    int half_size = round( size / 2.0);
    int i_size = half_size << 1;
    int top = round( int(y) - size/2);
    int left = round( int(x) - size/2);

    VW_ASSERT(left+i_size < (unsigned)integral.impl().cols(), 
	      vw::ArgumentErr() << "left out of bounds. "<< integral.impl().cols() <<" : "
	      << left+i_size << " [top left] " << top << " " << left << "\n");
    VW_ASSERT(top+i_size < (unsigned)integral.impl().rows(),
	      vw::ArgumentErr() << "top out of bounds. " << integral.impl().rows() <<" : "
	      << top+i_size << "\n");
    VW_ASSERT(left >= 0,
	      vw::ArgumentErr() << "left is to low. " << 0 << " : " << left << "\n");
    VW_ASSERT(top >= 0,
	      vw::ArgumentErr() << "top is to low. " << 0 << " : " << top << "\n");

    response = -integral.impl()(left, top);
    response += integral.impl()(left+i_size, top);
    response += 2*integral.impl()(left, top+half_size);
    response -= 2*integral.impl()(left+i_size, top+half_size);
    response -= integral.impl()(left, top+i_size);
    response += integral.impl()(left+i_size, top+i_size);
    
    return response;
  }

  // Horizontal Wavelet ( floating point arithmetic )
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate

  // Note: This Filter requires/recommends the use of an interpolated
  //       view of the integral
  template <class ViewT, class NumberT>
  typename boost::enable_if<boost::is_floating_point<NumberT>, float>::type
  inline HHaarWavelet( ImageViewBase<ViewT> const& integral,
			     NumberT const& x, NumberT const& y,
			     float const& size ) {

    VW_ASSERT( (integral.impl()(10,10) != integral.impl()(10.4,10)) ||
	       (integral.impl()(10,10) == integral.impl()(11,10) ),
	       vw::ArgumentErr() << "Input Integral doesn't seem to be interpolated" );
    VW_ASSERT( (integral.impl()(10,10) != integral.impl()(10,10.4)) ||
	       (integral.impl()(10,10) == integral.impl()(10,11) ),
	       vw::ArgumentErr() << "Input Integral doesn't seem to be interpolated" );

    float response;
    float half_size = size / 2.0;
    float top = float(y) - half_size;
    float left = float(x) - half_size;

    response = -integral.impl()(left, top);
    response += 2*integral.impl()(left+half_size, top);
    response -= integral.impl()(left+size, top);
    response += integral.impl()(left, top+size);
    response -= 2*integral.impl()(left+half_size, top+size);
    response += integral.impl()(left+size, top+size);
    
    return response;
  }

  // Vertical Wavelet ( floating point arithmetic )
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate

  // Note: This Filter requires/recommends the use of an interpolated
  //       view of the integral
  template <class ViewT, class NumberT>
  typename boost::enable_if<boost::is_floating_point<NumberT>, float>::type
  inline VHaarWavelet( ImageViewBase<ViewT> const& integral,
			     NumberT const& x, NumberT const& y,
			     float const& size ) {

    VW_ASSERT( (integral.impl()(10,10) != integral.impl()(10.4,10)) ||
	       (integral.impl()(10,10) == integral.impl()(11,10) ),
	       vw::ArgumentErr() << "Input Integral doesn't seem to be interpolated" );
    VW_ASSERT( (integral.impl()(10,10) != integral.impl()(10,10.4)) ||
	       (integral.impl()(10,10) == integral.impl()(10,11) ),
	       vw::ArgumentErr() << "Input Integral doesn't seem to be interpolated" );

    float response;
    float half_size = size / 2.0;
    float top = float(y) - half_size;
    float left = float(x) - half_size;

    response = -integral.impl()(left, top);
    response += integral.impl()(left+size, top);
    response += 2*integral.impl()(left, top+half_size);
    response -= 2*integral.impl()(left+size, top+half_size);
    response -= integral.impl()(left, top+size);
    response += integral.impl()(left+size, top+size);
    
    return response;
  }

}} // end namespace vw::ip

#endif // __VW_INTERESTPOINT_INTEGRALIMAGE_H__
