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


/// \file Descriptor.h
///
/// Basic classes and functions for generating interest point
/// descriptors with focus on using integral images. So this code
/// specifically does not make crops of the original image. Instead
/// descriptor are expected to rotate/scale their response themselves.
///
#ifndef __VW_INTERESTPOINT_INTEGRAL_DESCRIPTOR_H__
#define __VW_INTERESTPOINT_INTEGRAL_DESCRIPTOR_H__

#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/IntegralImage.h>

namespace vw {
namespace ip {

  template <class ImplT>
  class IntegralDescriptorGeneratorBase {

    ImageView<float> m_integral;

  public:

    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    // Given an image and a list of interest points, set the
    // descriptor field of the interest points using compute_descriptor()
    // method.
    template <class ViewT>
    void operator() ( ImageViewBase<ViewT> const& image, InterestPointList& points ) {
      (*this)( image, points.begin(), points.end() );
    }
    template <class ViewT, class IterT>
    void operator() ( ImageViewBase<ViewT> const& image,
                      IterT start, IterT end ) {

      // Timing
      Timer total("\tTotal elapsed time", DebugMessage, "interest_point");

      m_integral = IntegralImage(pixel_cast<float>(channel_cast<float>(image.impl())));

      for (IterT i = start; i != end; i++ ) {
        // Wrapping integral image for interpolation
        i->descriptor.set_size( impl().descriptor_size() );
        impl().compute_descriptor( interpolate(m_integral), *i );
      }
    }

    // Default suport size ( i.e. descriptor window)
    int support_size() { return 41; }
    // Default descriptor(vector) length
    int descriptor_size() { return 128; }
  };

  // Simple Scaled Gradient Descriptor (v2)
  struct SGrad2DescriptorGenerator : public IntegralDescriptorGeneratorBase<SGrad2DescriptorGenerator> {

    std::list<Vector3> m_box_properties;
    std::list<Vector2> m_histogram_samp;
    std::list<Vector2> m_33_samp;

    SGrad2DescriptorGenerator() {
      // Push back constant box properties
      m_box_properties.push_back(Vector3(-2,-2,1));
      m_box_properties.push_back(Vector3(2,-2,1));
      m_box_properties.push_back(Vector3(2,2,1));
      m_box_properties.push_back(Vector3(-2,2,1));
      m_box_properties.push_back(Vector3(-3.75,-3.75,1.25));
      m_box_properties.push_back(Vector3(0,-3.75,1.25));
      m_box_properties.push_back(Vector3(3.75,-3.75,1.25));
      m_box_properties.push_back(Vector3(3.75,0,1.25));
      m_box_properties.push_back(Vector3(3.75,3.75,1.25));
      m_box_properties.push_back(Vector3(0,3.75,1.25));
      m_box_properties.push_back(Vector3(-3.75,3.75,1.25));
      m_box_properties.push_back(Vector3(-3.75,0,1.25));
      m_box_properties.push_back(Vector3(-6.4,0,1.6));
      m_box_properties.push_back(Vector3(0,-6.4,1.6));
      m_box_properties.push_back(Vector3(6.4,0,1.6));
      m_box_properties.push_back(Vector3(0,6.4,1.6));
      // Push back histogram dest samp;
      m_histogram_samp.push_back(Vector2(0,1));
      m_histogram_samp.push_back(Vector2(0.707106781186547,0.707106781186547));
      m_histogram_samp.push_back(Vector2(1,0));
      m_histogram_samp.push_back(Vector2(0.707106781186547,-0.707106781186547));
      m_histogram_samp.push_back(Vector2(0,-1));
      m_histogram_samp.push_back(Vector2(-0.707106781186547,-0.707106781186547));
      m_histogram_samp.push_back(Vector2(-1,0));
      m_histogram_samp.push_back(Vector2(-0.707106781186547,0.707106781186547));
      // Push back 3x3 sample locations;
      m_33_samp.push_back(Vector2(-1,-1));
      m_33_samp.push_back(Vector2(-1, 0));
      m_33_samp.push_back(Vector2(-1, 1));
      m_33_samp.push_back(Vector2(0,-1));
      m_33_samp.push_back(Vector2(0, 0));
      m_33_samp.push_back(Vector2(0, 1));
      m_33_samp.push_back(Vector2(1,-1));
      m_33_samp.push_back(Vector2(1, 0));
      m_33_samp.push_back(Vector2(1, 1));
    }

    template <class ViewT>
    void compute_descriptor ( ImageViewBase<ViewT> const& integral,
                              InterestPoint& ip ) {
      uint32 write_idx = 0;
      float scale = ip.scale;
      float co = cos(ip.orientation);
      float si = sin(ip.orientation);
      Matrix2x2 rotate( co, -si, si, co );

      // Iterating through all sample locations for histograms
      for ( std::list<Vector3>::const_iterator center = m_box_properties.begin();
            center != m_box_properties.end(); center++ ) {

        uint32 s_write_idx = 0;
        std::vector<Vector2> responses(18);
        // Sampling 3x3 locations
        for ( std::list<Vector2>::const_iterator offset = m_33_samp.begin();
              offset != m_33_samp.end(); offset++ ) {
          Vector2 location = scale*rotate*(center->z()*(*offset)+subvector((*center),0,2));
          responses[s_write_idx] = rotate*Vector2( HHaarWavelet(integral.impl(),
                                                                location[0]+ip.x,
                                                                location[1]+ip.y,
                                                                2*scale*center->z()),
                                                   0);
          s_write_idx++;
          responses[s_write_idx] = rotate*Vector2( 0,
                                                   VHaarWavelet(integral.impl(),
                                                                location[0]+ip.x,
                                                                location[1]+ip.y,
                                                                2*scale*center->z()) );
          s_write_idx++;
        }

        // Fitting to histogram
        for ( std::list<Vector2>::const_iterator hist_bin = m_histogram_samp.begin();
              hist_bin != m_histogram_samp.end(); hist_bin++ ) {
          for ( size_t i = 0; i < responses.size(); i++ ) {
            double dot = dot_prod( *hist_bin, responses[i] );
            if ( dot > 0 )
              ip.descriptor[write_idx] += dot;
          }
          write_idx++;
        }
      }

      // Normalizing for lighting invariance
      float distance = norm_2( ip.descriptor );
      for ( size_t i = 0; i < 128; i++ )
        ip.descriptor[i] /= distance;

      // Threshold + Normalizing for non-linear lighting
      for ( size_t i = 0; i < 128; i++ )
        if ( ip.descriptor[i] > 0.2 )
          ip.descriptor[i] = 0.2;
      distance = norm_2( ip.descriptor );
      for ( size_t i = 0; i < 128; i++ )
        ip.descriptor[i] /= distance;
    }
  };

  // M-SURF Descriptor. This is an implementation of MU-SURF from
  // CenSurE's paper that includes the ability to take in consideration of
  // orientation. This is different from SURF in that samples are weighted
  // with a gaussian to reduce edge effects.

  // -- held for the moment --

}}

#endif//__VW_INTERESTPOINT_INTEGRAL_DESCRIPTOR_H__
