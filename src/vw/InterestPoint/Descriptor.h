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

/// \file Descriptor.h
/// 
/// Basic classes and functions for generating interest point descriptors.
/// 
#ifndef __INTERESTPOINT_DESCRIPTOR_H__
#define __INTERESTPOINT_DESCRIPTOR_H__

#include <vector>
#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Image/Transform.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/InterestPoint/ImageOctave.h>
#include <vw/InterestPoint/InterestData.h>

namespace vw { 
namespace ip {

  /// This metric can be used to measure the error between a keypoint
  /// p2 and a second keypoint p1 that is transformed by a 3x3 matrix
  /// H.  This is predominately used when matching keypoints using
  /// RANSAC.
  struct KeypointErrorMetric {
    template <class RelationT, class ContainerT>
    double operator() (RelationT const& H,
                       ContainerT const& p1, 
                       ContainerT const& p2) const {
      return vw::math::norm_2( Vector3(p2.x,p2.y,1) - H * Vector3(p1.x,p1.y,1));
    }
  };


  /// Select only the interest points that fall within the specified bounding box.
  template <class RealT>
  KeypointList crop(KeypointList const& interest_points, BBox<RealT,2> const& bbox) {
    KeypointList return_val;
    for (KeypointList::iterator i = interest_points.begin(); i != interest_points.end(); ++i) {
      if (bbox.contains(Vector<RealT,2>(RealT((*i).x), RealT((*i).y))))
        return_val.push_back(*i);
    }
    return return_val;
  }

#define SUPPORT_SIZE 41
  
  /// Get the size x size support region around an interest point,
  /// rescaled by the scale factor and rotated by the specified
  /// angle.
  /// transform takes care of interpolation and edge extension.
  template <class ViewT>
  inline ImageView<typename ViewT::pixel_type>
  get_support( float x, float y, float scale, float ori,
               ImageViewBase<ViewT> const& source, int size=SUPPORT_SIZE ) {
    float half_size = ((float)(size - 1)) / 2.0f;
    float scaling = 1.0f / scale;
    // This is mystifying - why won't the four-arg compose work?
    return transform(source.impl(),
                     compose(TranslateTransform(half_size, half_size),
                     compose(ResampleTransform(scaling, scaling),
                             RotateTransform(-ori),
                             TranslateTransform(-x, -y))),
                     size, size);
  }

  /// Get the support region around an interest point, scaled and
  /// rotated appropriately.
  template <class ViewT>
  inline ImageView<typename ViewT::pixel_type>
  get_support( InterestPoint const& pt, ImageViewBase<ViewT> const& source,
               int size=SUPPORT_SIZE ) {
    return get_support(pt.x, pt.y, pt.scale, pt.orientation, source.impl(), size);
  }

//   /// CRTP base class for descriptor generating methods.
//   /// A descriptor generator subclass can easily be created by
//   /// implementing cache_support and compute_descriptor_from_support,
//   /// or compute_descriptors can be reimplemented entirely.
//   template <class ImplT>
//   struct DescriptorBase {
//     /// Returns the derived implementation type.
//     ImplT& impl() { return *static_cast<ImplT*>(this); }

//     /// Returns the derived implementation type.
//     ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

//     /// Compute descriptors for a set of interest points from the given
//     /// source (an image or set of images).
//     inline void compute_descriptors( KeypointList& points ) const {
//       for (KeypointList::iterator i = points.begin(); i != points.end(); ++i) {
//         impl()(*i);
//       }
//     }
//   };


  template <class ViewT, class DescriptorT> 
  void compute_descriptors( ImageViewBase<ViewT> const& image, KeypointList& points, DescriptorT const& descriptor ) {
    for (KeypointList::iterator i = points.begin(); i != points.end(); ++i) {
      // First we ompute the support region based on the keypoint 
      ImageView<typename ViewT::pixel_type> support = get_support(*i, image.impl());

      // Next, we pass the support region and the keypoint to the
      // descriptor generator.
      i->descriptor = descriptor(support);
    }
  }


  /// A basic example descriptor class. The descriptor for an interest
  /// point is the support region around the point. It is normalized
  /// to provide some tolerance to changes in illumination.
  class PatchDescriptor {
  public:
    template <class ViewT>
    Vector<float> operator() (ImageViewBase<ViewT> const& support) const {
      Vector<float> result;
      result.set_size(SUPPORT_SIZE * SUPPORT_SIZE);
      for (int j = 0; j < SUPPORT_SIZE; j++)
        for (int i = 0; i < SUPPORT_SIZE; i++)
          result(i*SUPPORT_SIZE + j) = support.impl()(i,j);
      return normalize(result);
    }
  };

  // TODO: assign_patch_descriptors convenience

}} // namespace vw::ip

#endif //__INTERESTPOINT_DESCRIPTOR_H__
