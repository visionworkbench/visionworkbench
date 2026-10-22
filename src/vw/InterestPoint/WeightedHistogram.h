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


/// \file WeightedHistogram.h
///
/// A helper for constructing a 2D Gaussian kernel.
///
#ifndef _WEIGHTED_HISTOGRAM_H_
#define _WEIGHTED_HISTOGRAM_H_

#include <vw/Image/ImageView.h>
#include <vw/Image/Filter.h>
#include <vector>

namespace vw { namespace ip {

/// Construct a 2D Gaussian kernel. This is templatized so that kernel
/// could be any 2D array that provides the appropriate accessors and
/// resize method.  Note that for Gaussian filtering of images, the
/// separable method is preferred with two 1D kernels.
/// The kernel size may be specified; if it is not then the default
/// value of zero indicates the kernel size should be computed from sigma.
template <class KernelT>
int make_gaussian_kernel_2d( KernelT& kernel, float sigma, int usewidth=0 ) {
  std::vector<float> kernel_1d;
  generate_gaussian_kernel(kernel_1d, sigma, usewidth);
  int kerwidth = kernel_1d.size();
  kernel.set_size(kerwidth,kerwidth);

  // Put in Gaussian values
  for (int j=0; j<kerwidth; j++)
    for (int i=0; i<kerwidth; i++)
      kernel(i,j) = kernel_1d[i] * kernel_1d[j];

  return 0;
}

}} // namespace vw::ip

#endif
