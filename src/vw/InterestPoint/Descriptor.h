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
  std::vector<InterestPoint> crop(std::vector<InterestPoint> const& interest_points, BBox<RealT,2> const& bbox) {
    std::vector<InterestPoint> return_val;
    for (int i = 0; i < interest_points.size(); ++i) {
      if (bbox.contains(Vector<RealT,2>(RealT(interest_points[i].x), RealT(interest_points[i].y))))
        return_val.push_back(interest_points[i]);
    }
    return return_val;
  }

#define SUPPORT_SIZE 41
  
  /// Get the size x size support region around an interest point,
  /// rescaled by the scale factor and rotated by the specified
  /// angle.
  /// transform takes care of interpolation and edge extension.
  template <class T>
  int inline get_support( ImageView<T>& support,
		   float x, float y, float scale, float ori,
		   const ImageView<T>& source, int size=SUPPORT_SIZE ) {
    float half_size = ((float)(size - 1)) / 2.0f;
    float scaling = 1.0f / scale;
    // This is mystifying - why won't the four-arg compose work?
    support = transform(source,
			compose(TranslateTransform(half_size, half_size),
			compose(ResampleTransform(scaling, scaling),
				RotateTransform(-ori),
				TranslateTransform(-x, -y))),
			size, size);

    return 0;
  }

  /// Get the support region around an interest point, scaled and
  /// rotated appropriately.
  template <class T>
  int inline get_support( ImageView<T>& support,
			  const InterestPoint& pt,
			  const ImageView<T>& source,
			  int size=SUPPORT_SIZE ) {
    return get_support(support, pt.x, pt.y, pt.scale, pt.orientation,
		       source, size);
  }

  /// CRTP base class for descriptor generating methods.
  /// A descriptor generator subclass can easily be created by
  /// implementing cache_support and compute_descriptor_from_support,
  /// or compute_descriptors can be reimplemented entirely.
  template <class ImplT, class SourceT>
  struct DescriptorBase {
    /// Returns the derived implementation type.
    ImplT& impl() { return *static_cast<ImplT*>(this); }

    /// Returns the derived implementation type.
    ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

    /// Compute descriptors for a set of interest points from the given
    /// source (an image or set of images).
    int compute_descriptors( std::vector<InterestPoint>& points,
			    const SourceT& source ) {
      for (int i = 0; i < points.size(); i++) {
        impl().cache_support(points[i], source );
        impl().compute_descriptor_from_support(points[i]);
      }
      return 0;
    }
  };

  /// Generate descriptors for the set of interest points using the
  /// given source data and descriptor generator.
  template <class DescriptorT, class SourceT>
  int generate_descriptors(std::vector<InterestPoint>& points,
			   const SourceT& source,
			   DescriptorBase<DescriptorT, SourceT>& desc_gen) {
    return desc_gen.impl().compute_descriptors(points, source);
  }

  /// A basic example descriptor class. The descriptor for an interest
  /// point is the support region around the point. It is normalized
  /// to provide some tolerance to changes in illumination.
  template <class T>
  class PatchDescriptor : public DescriptorBase<PatchDescriptor<T>, ImageView<T> >
  {
    ImageView<T> support;

  public:
    typedef T real_type;

    inline int cache_support(InterestPoint& pt,
		             const ImageView<T>& source) {
      return vw::ip::get_support(support, pt, source); // 41x41
    }

    int compute_descriptor_from_support(InterestPoint& pt) {
      pt.descriptor.set_size(SUPPORT_SIZE * SUPPORT_SIZE);
      for (int i = 0; i < SUPPORT_SIZE; i++)
	for (int j = 0; j < SUPPORT_SIZE; j++)
	  pt.descriptor(i*SUPPORT_SIZE + j) = support(i,j);
      pt.descriptor = normalize(pt.descriptor);

      return 0;
    }
  };

}} // namespace vw::ip

#endif //__INTERESTPOINT_DESCRIPTOR_H__
