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

/// \file WeightedHistogram.h
/// 
/// Functions for constructing histograms from pixels in a given image.
/// 
#ifndef _WEIGHTED_HISTOGRAM_H_
#define _WEIGHTED_HISTOGRAM_H_

#include <vw/Image/ImageView.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Manipulation.h>
#include <vector>
#include <algorithm>
#include <math.h>

// Yet another magic number from Lowe.
#define INTEREST_POINT_MODE_THRESHOLD (0.8)

namespace vw { namespace ip {

/// Construct histogram from pixels in given image
template <class T1>
inline
void weighted_histogram(const ImageView<T1>& val_image,
                        const ImageView<T1>& weight_image,
                        std::vector<double>& histo,
                        double min, double max, unsigned n_bins) {
  // Make sure images are same size
  assert( val_image.cols() == weight_image.cols() );
  assert( val_image.rows() == weight_image.rows() );

  histo.resize(n_bins);
  std::fill(histo.begin(),histo.end(),0.0);
  double x0 = double(min);
  double s = double(n_bins-1)/(double(max)-x0);

  typename ImageView<T1>::iterator weight_iter = weight_image.begin();
  for (typename ImageView<T1>::iterator iter = val_image.begin(); iter != val_image.end(); ++iter, ++weight_iter) {
    int index = int(0.5+s*(double(*iter) - x0));
    
    if (index>=0 && (unsigned)index<n_bins) histo[index] += *weight_iter;
  }
}

/// Construct a 2D Gaussian kernel. This is templatized so that kernel
/// could be any 2D array that provides the appropriate accessors and
/// resize method.  Note that for Gaussian filtering of images, the
/// separable method is preferred with two 1D kernels.
/// The kernel size may be specified; if it is not then the default
/// value of zero indicates the kernel size should be computed from
/// sigma.
template <class KernelT>
int make_gaussian_kernel_2d( KernelT& kernel, float sigma, int usewidth=0 ) {
  std::vector<float> kernel_1d;
  generate_gaussian_kernel(kernel_1d, (double)sigma, usewidth);
  int kerwidth = kernel_1d.size();
  kernel.set_size(kerwidth,kerwidth);
  int halfwidth = (kerwidth-1)/2;
  
  // Put in Gaussian values
  for (int j=0; j<kerwidth; j++)
    for (int i=0; i<kerwidth; i++)
      kernel(i,j) = kernel_1d[i] * kernel_1d[j];

  return 0;
}

/// Compute (gaussian weight) * (edge magnitude) kernel.
template <class T>
void weighted_magnitude(ImageView<T>& weight, const ImageView<T>&
                        mag_image, int x, int y, T sigma,
                        int usewidth=0) {
  ImageView<T> kernel;
  make_gaussian_kernel_2d(kernel, sigma, usewidth);
  int width = kernel.rows();
  int halfwidth = (width - 1) / 2;

  weight = crop(mag_image, x - halfwidth, y - halfwidth, width,
                width) * kernel;
}

/// Adds a sample to the orientation histogram using linear
/// interpolation.
template <class T>
void orientation_interpolate(std::vector<T>& histo,
                             unsigned n_bins, T ori, T mag) {
  T bin_size = (T)2 * M_PI / (T)n_bins;
  T ratio = (ori - bin_size) / bin_size;
  int n = (int)floor((ori + M_PI) / bin_size);
  histo[n] += ((T)1 - ratio) * mag;
  histo[(n+1) % n_bins] += ratio * mag;
}

/// Construct orientation histogram for a region
template <class T>
void orientation_histogram(const ImageView<T>& ori_image,
                           const ImageView<T>& mag_image,
                           std::vector<T>& histo, int x, int y,
                           T sigma, unsigned n_bins) {
  ImageView<T> weight;
  weighted_magnitude(weight, mag_image, x, y, sigma);
  int width = weight.rows();
  int halfwidth = (width - 1) / 2;
  histo.resize(n_bins);
  std::fill(histo.begin(),histo.end(),0.0);

  for (int j = x - halfwidth; j <= x + halfwidth; j++)
    for (int i = y - halfwidth; i <= y + halfwidth; i++)
      orientation_interpolate(histo, n_bins, ori_image(x,y),
                              mag_image(x,y));
}

/// Suppresses non-maxima in a vector
template <class T>
void non_max_suppression( std::vector<T>& hist,
			  bool wrap=false ) {
  // Copy input to a buffer so we can check values in one while
  // modifying the other
  std::vector<T> buf(hist);
  int numel = hist.size();

  // Left edge
  if (hist[0]<hist[1]) buf[0] = 0;
  // Center
  for (int i=1; i<numel-1; ++i){
    if (hist[i]<hist[i+1]) buf[i] = 0;
    if (hist[i]<hist[i-1]) buf[i] = 0;
  }
  // Right edge
  if (hist[numel-1]<hist[numel-2]) buf[numel-1] = 0;

  // Finally, if we are wrapping, check left edge vs. right edge
  if (wrap){
    if (hist[0]<hist[numel-1]) buf[0] = 0;
    if (hist[0]>hist[numel-1]) buf[numel-1] = 0;
  }

  // Copy the non-max suppressed buffer back into hist
  hist = buf;

}

/// Correlate a 1D vector with a 1D kernel, wrapping the signal.  This
/// is useful for kernel density estimates on histograms of angles,
/// where the first and last histogram bin are adjacent in angle space.
template <class T1, class T2>
int filter_1d( std::vector<T1>& vec,
	       std::vector<T2>& kernel,
	       bool wrap=false ) {
  int numel = vec.size();
  int kerwidth = kernel.size();
  assert( kerwidth%2 ); // make sure kernel width is odd
  int halfwidth = (kerwidth-1)/2;

  // Temp buffer large enough to hold one col of the image
  std::vector<T1> buf(numel);
  fill(buf.begin(), buf.end(), 0.0);
  
  // Correlate each element with the kernel
  for (int i=0; i<numel; i++){
    // correlate elements with the kernel centered on element i
    for (int di=-halfwidth; di<=halfwidth; di++){
      // kernel index:
      int ki = halfwidth+di;
      // vector index:
      if (wrap){
	// If we are wrapping, then wrap the vector index
	int vi = (i+di+numel)%numel;
	//printf( "wrap: i= %d vi= %d ki= %d\n", i, vi, ki );
	buf[i] = buf[i] + vec[vi] * kernel[ki];
      } else {
	// If we are not wrapping, truncate the vector index
	int vi = i+di;
	if ((vi>=0)&&(vi<numel))
	  buf[i] = buf[i] + vec[vi] * kernel[ki];
      }
    }
  }
  // Copy buffer into the vector
  for (int i=0; i<numel; i++){
    vec[i] = buf[i];
  }

  return 0;
}

/// Smooth a histogram with a Gaussian kernel
template <class T>
int smooth_weighted_histogram( std::vector<T>& histo,
			       float sigma) {
  std::vector<float> kernel;
  vw::generate_gaussian_kernel( kernel, sigma );
  bool wrap = true;
  filter_1d( histo, kernel, wrap );
  
  return 0;
}

/// Finds the max bin in a histogram.  Mode is a vector in case there
/// is a second equal or nearly equal bin in which case both are returned.
template <class T>
void find_weighted_histogram_mode( const std::vector<T>& hist_in,
                                   std::vector<int>& mode) {
  mode.clear();
  
  // Make a copy so we can destroy it in non_max_suppression as we go
  // along
  std::vector<T> hist = hist_in;
  
  // Suppress points other than local maxima
  bool wrap=true;
  non_max_suppression( hist, wrap );
  
  // Find max value & location.
  double max_val=0.0;
  int max_bin=0;
  for (unsigned int i=0; i<hist.size(); i++){
    if (hist[i]>max_val){
      max_val = hist[i];
      max_bin = i;
    }
  }
  
  // Put max into the mode vector
  mode.push_back( max_bin );
  
  // Zero the max.
  hist[max_bin] = 0;

  // Find other local maxima (modes) within a specified threshold of
  // the global maximum.  Don't include too many...
  while ( ( max_val > INTEREST_POINT_MODE_THRESHOLD*mode[0] ) &&
	  ( mode.size()<5 ) ){
    max_val=0.0;
    max_bin=0;
    for (unsigned int i=0; i<hist.size(); i++){
      if (hist[i]>max_val){
        max_val = hist[i];
        max_bin = i;
      }
    }
    // If next largest bin is within a threshold of the first, add it
    // to the vector of modes
    if ( max_val > INTEREST_POINT_MODE_THRESHOLD*mode[0] ) mode.push_back( max_bin );
  }
  
}

}} // namespace vw::ip

#endif 

