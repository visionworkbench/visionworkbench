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

/// \file SIFT_Descriptor.h
/// 
/// Class for generating SIFT descriptors.
/// 
#ifndef __INTERESTPOINT_SIFTDESCRIPTOR_H__
#define __INTERESTPOINT_SIFTDESCRIPTOR_H__

#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/WeightedHistogram.h>

#include <stdio.h>

namespace vw { 
namespace ip {

// Magic number defined by Lowe.
#define SIFT_MAGNITUDE_THRESHOLD 0.2

/// A SIFT descriptor is formed by taking a 16x16 support region
/// around the interest point, dividing it into a 4x4 grid, and
/// creating an 8-bin orientation histogram for each cell of the
/// grid. The descriptor is a vector of length 4x4x8 = 128.
///
/// The source data is the ImageOctaveHistory, which can be
/// recorded using ScaledInterestPointDetector::record_history.
  template <class T>
  class SIFT_Descriptor : public DescriptorBase<SIFT_Descriptor<T>, ImageOctaveHistory<ImageInterestData<T> > >
  {
    ImageView<T> mag_region;
    ImageView<T> ori_region;
    ImageView<T> kernel;
    T bin_size, two_pi;

  public:
    typedef T real_type;
    typedef ImageOctaveHistory<ImageInterestData<T> > source_type;

    SIFT_Descriptor() {
      make_gaussian_kernel_2d(kernel, 9.0, 16);
      bin_size = M_PI / (T)4;
      two_pi = M_PI * (T)2;
    }

    inline int cache_support(InterestPoint& pt,
                 const source_type& source) {
      // Select level of Gaussian blur from scale.
      ImageInterestData<T> data = source.image_at_scale(pt.scale);
      // Get support regions, without performing any further rescaling.
      get_support(mag_region, pt.x, pt.y, 1.0, pt.orientation, data.mag, 16);
      get_support(ori_region, pt.x, pt.y, 1.0, pt.orientation, data.ori, 16);
      // Rotate gradient orientations by point orientation
      ori_region += (-pt.orientation);
    }

    int compute_descriptor_from_support(InterestPoint& pt) {
      pt.descriptor.set_size(128);

      // Apply gaussian weighting to magnitudes
      ImageView<T> weight_region = mag_region * kernel;

      for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
          // Get orientation between 0 and 2*pi
          T ori = ori_region(i, j);
          while (ori < 0) ori += two_pi;
          while (ori >= two_pi) ori -= two_pi;

          T weight = weight_region(i, j);
          int bin = (int)(ori / bin_size);

          // Tri-linear interpolation
          int a = (i - 2) / 4;
          int b = (j - 2) / 4;
          // Weight will be split across four histograms:
          // H(a,b), H(a+1,b), H(a,b+1), and H(a+1,b+1)
          for (int c = a; c < a + 2; c++) {
            for (int d = b; d < b + 2; d++) {
              if ((c >= 0) && (c < 4) && (d >= 0) && (d < 4)) {
                T weight_slice = (1.0 - fabs(i - 4.0*c - 1.5) / 4.0) *
                                 (1.0 - fabs(j - 4.0*d - 1.5) / 4.0);
                // Interpolate between closest two orientation bins
                T r = (ori - (T)bin * bin_size) / bin_size;
                int hist_offset = 32 * d + 8 * c;
                pt.descriptor(hist_offset + bin) += (1.0 - r) * weight_slice;
                pt.descriptor(hist_offset + (bin+1)%8) += r * weight_slice;
              }
            }
          }
        }
      }

      // Normalize to reduce the effects of illumination change
      pt.descriptor = normalize(pt.descriptor);

      // Reduce influence of large gradient magnitudes
      for (int i = 0; i < pt.descriptor.size(); i++) {
        if (pt.descriptor(i) > SIFT_MAGNITUDE_THRESHOLD)
          pt.descriptor(i) = SIFT_MAGNITUDE_THRESHOLD;
        if (pt.descriptor(i) < 0)
          pt.descriptor(i) = 0;
      }
      pt.descriptor = normalize(pt.descriptor);

      return 0;
    }
  };

}} // namespace vw::ip

#endif //__INTERESTPOINT_SIFTDESCRIPTOR_H__
