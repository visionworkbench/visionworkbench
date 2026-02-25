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


/// \file ImageOctave.h
///
/// ImageOctave class for constructing the scale space of a source image.
///
#ifndef __VW_INTERESTPOINT_IMAGEOCTAVE_H__
#define __VW_INTERESTPOINT_IMAGEOCTAVE_H__

// Vision Workbench
#include <vw/Image/ImageView.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageIO.h>

// STL
#include <vector>
#include <iostream>

// Boost
#include <boost/static_assert.hpp>

// Magic numbers. Learn these parameters?
#define INITIAL_SIGMA (1.6)
#define CAMERA_SIGMA (0.5)

namespace vw { namespace ip {

/// This class provides functionality to construct an octave of images
/// from one image, or to recursively (and destructively - use also
/// ImageOctaveHistory to keep the entire pyramid) construct the
/// next octave from an existing octave.  This octave makes up one part
/// of a scale space pyramid.  Specifically, it makes up the part of
/// scale space that covers one doubling in scale, where the number of
/// planes is dictated by the scale ratio between planes (or vice
/// versa).  In addition, since we want to bound peaks in all three
/// (x,y,s) directions, it constructs one plane above and one plane
/// below the octave. If the scale ratio is 2^(1/P), the octave
/// contains P+2 planes.
template <class ViewT = ImageView<float> >
class ImageOctave {
  BOOST_STATIC_ASSERT(IsImageView<ViewT>::value);

public:
  typedef ViewT scale_type;

  int base_scale;            // initially = 1, then doubles every build_next()
  int scales_per_octave;     // P, where scale ratio = 2^(1/P)
  int num_planes;            // Number of planes (P+2) for ratio of 2^(1/P)
  float init_sigma;          // Initial sigma to apply to base level image
  float sigma_ratio;         // Ratio of sigmas between levels:
  std::vector<float> sigma;  // Sigmas corresponding to scales
  std::vector<ViewT> scales; // Scaled images in the octave

  /// This constructor is intended for building the first octave from a
  /// source image.
  template <class ViewT_in>
  ImageOctave( ImageViewBase<ViewT_in> const& src_im, int numscales) {
    base_scale = 1;
    set_scales( numscales );

    build_using( src_im );
  }

  /// Calculate the scale corresponding to a plane index.
  float plane_index_to_scale( float plane_index ) const {
    float s = base_scale*pow(2.0f, plane_index / ((float)scales_per_octave) );
    return s;
  }

  /// Static function for calculating plane index from scale.
  static int scale_to_plane_index( int base, int scales, float scale) {
    return (int)(scales * (log(scale) - log((float)base))/M_LN2 + 0.00001);
  }

  /// Calculate the plane index most closely matching a scale.
  int scale_to_plane_index( float scale ) const {
    return ImageOctave::scale_to_plane_index(base_scale, scales_per_octave, scale);
  }

  /// Sets the number of scales, the initial sigma, the ratio of
  /// sigmas, and the vector of sigmas.
  int set_scales( int scalesperoctave ) {
    // Set number of scales per octave
    scales_per_octave = scalesperoctave;

    // Set number of planes, which adds one plane above and below
    num_planes = scales_per_octave + 2;

    // Initial sigma
    init_sigma = INITIAL_SIGMA;

    // Compute sigma ratio
    sigma_ratio = pow(2.0,1.0/((float)scales_per_octave));

    // Set up sigmas
    sigma.resize(num_planes);
    sigma[0] = init_sigma / sigma_ratio;
    for (int i=1; i<num_planes; i++){
      sigma[i] = sigma[i-1] * sigma_ratio;
    }
    return 0;
  }

  /// Build image octave using designated source image.  The number of
  /// scales and the sigmas should already be computed.
  // TODO: Use supersampled view as the first scale
  // TODO: this could be done a little more efficiently
  template <class ViewT_in>
  int build_using(ImageViewBase<ViewT_in> const& src_im ) {
    // Construct enough image plane containers for the number of
    // scales in the octave
    scales.clear();
    scales.reserve(num_planes);

    // Scales constructed by blurring.  Assume some sigma for the first one.
    vw::vw_out(VerboseDebugMessage, "interest_point") << "\tAssuming camera_sigma " << CAMERA_SIGMA
                                    << std::endl;
    float camera_sigma = CAMERA_SIGMA;

    // Sigma to use for blurring step to achieve a final sigma in each
    // plane.  With repeated Gaussian blurring, the sigmas add in
    // quadrature, so a few smaller kernels in succession can produce
    // the output of a larger kernel with fewer operations.
    if (sigma[0]>camera_sigma){
      float use_sigma = sqrt( sigma[0]*sigma[0] - camera_sigma*camera_sigma );
      vw::vw_out(VerboseDebugMessage, "interest_point") << "\tMaking plane " << 0 << " using sigma "
                                      << use_sigma << " so final sigma is "
                                      << sigma[0] << std::endl;
      scales.push_back(ViewT(vw::gaussian_filter(src_im, use_sigma)));
    } else {
      scales.push_back(ViewT(src_im.impl()));
    }

    // Each next plane is blurred version of the previous
    for (int i=1; i<num_planes; i++){
      float use_sigma = sqrt( sigma[i]*sigma[i] - sigma[i-1]*sigma[i-1] );
      vw::vw_out(VerboseDebugMessage, "interest_point") << "\tMaking plane " << i << " using sigma "
                                      << use_sigma << " so final sigma is "
                                      << sigma[i] << std::endl;
      scales.push_back(ViewT(vw::gaussian_filter(scales[i-1], use_sigma)));
    }

    return 0;
  }

  /// Build the next octave using the current one.  The number of
  /// scales and the sigmas should already be computed.
  int build_next() {
    // Keep track of base scale.  For the first octave, it is 1.  Each
    // build_next() call multiplies it by two.  This number multiplied
    // by the sigma for each plane is the scale relative to the
    // original image.  The sigma for each plane is the scale relative
    // to the base image in the current octave.
    base_scale *= 2;

    // Plane 0 is the current plane #(N-2), downsampled by 2.  For
    // example, if there are 3 planes per octave, and the two extra
    // planes for bounding points in scale space, then the
    // source->dest images go like:
    // old octave:  0 1 2 3 4
    // new octave:        0 1 2 3 4
    // so that image #3 in the first octave is subsampled to produce
    // #0 in the new octave, image #4 is subsampled to produce #1 in
    // the new octave, and the rest are built by successive blurring
    // operations.
    int source_for_new_0 = num_planes-2;
    scales[0] = ViewT(vw::subsample(scales[source_for_new_0], 2));

    // Plane 1 is also the subsampled plane #(N-1)
    int source_for_new_1 = num_planes-1;
    scales[1] = ViewT(vw::subsample(scales[source_for_new_1], 2));

    // Now the rest (#2 and up) are blurred versions as before
    float use_sigma;
    for (int k=2; k<num_planes; k++){
      use_sigma = sqrt( sigma[k]*sigma[k] - sigma[k-1]*sigma[k-1] );
      vw::vw_out(VerboseDebugMessage, "interest_point") << "\tMaking plane " << k << " using sigma "
                                      << use_sigma << " so final sigma is "
                                      << sigma[k] << std::endl;
      scales[k] = ViewT(vw::gaussian_filter(scales[k-1], use_sigma));
    }

    return 0;
  }


  /// Save image octave to a set of image files.
  void write_images() const
  {
    char fname[256]; // output filename
    int imagenum; // which image number for unique filenames

    // Save each image plane in the octave
    for (int k=0; k<num_planes; k++){
      imagenum = (int)(log((float)base_scale)/log(2.0)) * num_planes + k;
      snprintf(fname, sizeof(fname), "scale_%02d.jpg", imagenum);
      vw::write_image( fname, scales[k] );
    }
  }
};

}} // namespace vw::ip

#endif   // _IMAGE_OCTAVE_H_
