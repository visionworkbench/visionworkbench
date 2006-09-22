// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file Filter.tcc
/// 
/// Two-dimensional image filter classes and functions.
/// 
/// These are the definitions of the non-trivial functions 
/// and methods declared in Filter.h.  See that file for 
/// more discussion.
///
#ifndef __VW_IMAGE_FILTER_TCC__
#define __VW_IMAGE_FILTER_TCC__

#include <vw/Core/Exception.h>
#include <vw/Math/Functions.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Filter.h>

/// Compute a Gaussian kernel.  The default size is seven times sigma 
/// rounded down to the nearest odd integer, or 3, whichever is larger.
template <class KernelT>
void vw::generate_gaussian_kernel( std::vector<KernelT>& kernel, double sigma, int size )
{
  if( sigma == 0 ) {
    kernel.clear();
    return;
  }
  if( size == 0 ) {
    size = (int)(7*sigma);
    if( size<3 ) size = 3;
    else if( size%2==0 ) size -= 1;
  }
  kernel.resize( size );

  unsigned center = size / 2; 
  double sum = 0.0, tap; 
  const double z = 1 / (sqrt(2.0) * sigma);

  // Even length filter.  Off center.
  if (size % 2 == 0) {
    for (unsigned i=0 ; i < center ; ++i) {
      tap = vw::erf((i+1.0) * z) - vw::erf(i * z);
      sum += tap;
      kernel[center+i] = kernel[center-i-1] = (KernelT)tap;
    }
    sum *= 2.0;
  }

  // Odd length filter.  Perfectly centered.
  else {
    for (unsigned i=1 ; i<=center; ++i) {
      tap = vw::erf((i+0.5) * z) - vw::erf((i-0.5) * z);
      sum += tap;
      kernel[center+i] = kernel[center-i] = (KernelT)tap;
    }
    sum *= 2.0;
    tap = vw::erf(0.5 * z) - vw::erf(-0.5 * z);
    sum += tap;
    kernel[center] = (KernelT)tap;
  }
	
  // Normalize the result
  assert(sum >= 0.0);
  double norm = 1.0 / sum;
  std::transform(kernel.begin(), kernel.end(), kernel.begin(),
                 std::bind2nd(std::multiplies<KernelT>(), (KernelT)norm));
}


// Compute a differentiation kernel.
//
// Assume that only the n lowest-order terms of the Taylor expansion
// of the function are present, where n is the size of the kernel.
// Build up a Taylor expansion matrix T(i,j)=(i-halfsize)^j/j!.  The
// rows of the inverse of this matrix are the first n differentation
// operators; just pick the one you want.  This should get reworked 
// if we ever want to properly support integer arithmetic.
template <class KernelT>
void vw::generate_derivative_kernel( std::vector<KernelT>& kernel, int deriv, int size )
{
  // Disable filter for zeroth derivative
  if( deriv == 0 ) {
    kernel.clear();
    return;
  }

  // Check and configure kernel size
  int minsize = deriv + (deriv%2) + 1;
  if( size == 0 ) size = minsize;
  else if( size < minsize ) throw ArgumentErr( "Derivative kernel too small for requested differentiation operator!" );
  else if( size%2 == 0 ) throw ArgumentErr( "Kernel must have odd dimensions!" );
  kernel.resize( size );

  // Test for most common cases
  if( deriv==1 && size==3 ) {
    kernel[0] = 0.5;
    kernel[1] = 0.0;
    kernel[2] = -0.5;
    return;
  }
  if( deriv==2 && size==3 ) {
    kernel[0] = 1.0;
    kernel[1] = -2.0;
    kernel[2] = 1.0;
    return;
  }

  // Compute the size-th order polynomial matrix
  Matrix<KernelT> pmat(size,size);
  int half_size = size/2;
  for( int j=0; j<size; ++j ) {
    int x = half_size - (int)j;
    KernelT term = 1;
    for( int i=0; i<size; ++i ) {
      pmat(i,j) = term;
      term *= x;
      term /= i+1;
    }
  }

  Vector<KernelT> dsel(size);
  dsel[deriv] = 1.0;
  Vector<KernelT> kv = inverse(pmat)*dsel;
  for( int i=0; i<size; ++i )
    kernel[i] = kv[i];
}


// Compute a two-dimensional Gaussian derivative kernel.
template <class KernelT>
void vw::generate_gaussian_derivative_kernel( ImageView<KernelT>& kernel, double sigma1, int deriv1, double sigma2, int deriv2, double angle, int size ) {
  kernel.set_size( size, size, 1 );
  double ca=vw::cos(angle), sa=vw::sin(angle), half=size/2;
  double scalar = 2*M_PI*sigma1*sigma2*vw::pow(-sigma1*sigma1,deriv1)*vw::pow(-sigma2*sigma2,deriv2);
  double sum = 0;
  for( int i=0; i<size; ++i ) {
    for( int j=0; j<size; ++j ) {
      double x=ca*(i-half)+sa*(j-half);
      double y=-sa*(i-half)+ca*(j-half);
      kernel(i,j) = vw::exp(-x*x/(2*sigma1*sigma1)) * vw::exp(-y*y/(2*sigma2*sigma2)) / scalar;
      if( deriv1==1 ) kernel(i,j) *= x;
      else if( deriv1==2 ) kernel(i,j) *= (x*x-sigma1*sigma1);
      if( deriv2==1 ) kernel(i,j) *= y;
      else if( deriv2==2 ) kernel(i,j) *= (y*y-sigma2*sigma2);
      sum += kernel(i,j);
    }
  }
  if( deriv1==0 && deriv2==0 ) kernel /= sum;
  else kernel -= sum/(size*size);
}


// Compute a two-dimensional Laplacian of Gaussian kernel.
template <class KernelT>
void vw::generate_laplacian_of_gaussian_kernel( ImageView<KernelT>& kernel, double sigma, int size ) {
  kernel.set_size( size, size, 1 );
  double half=size/2;
  double scalar = 2*M_PI*sigma*sigma*sigma*sigma*sigma*sigma;
  double sum = 0;
  for( int i=0; i<size; ++i ) {
    for( int j=0; j<size; ++j ) {
      double x=i-half, y=j-half;
      kernel(i,j) = vw::exp(-(x*x+y*y)/(2*sigma*sigma)) * (x*x+y*y-2*sigma*sigma) / scalar;
      sum += kernel(i,j);
    }
  }
  kernel -= sum/(size*size);
}

#endif  // __VW_IMAGE_FILTER_TCC__
