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

/// \file Image.h
/// 
/// Basic classes and structures for storing image interest points.
/// 
#ifndef __INTERESTPOINT_DESCRIPTOR_H__
#define __INTERESTPOINT_DESCRIPTOR_H__

#include <vector>
#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Image/Transform.h>
#include <vw/InterestPoint/ImageOctave.h>

namespace vw { 
namespace ip {

  /// A class for storing information about an interest point.
  struct InterestPoint {

    /// subpixel (col,row) location of point
    float x,y;

    /// scale of point.  This may come from the pyramid level, from
    /// interpolating the interest function between levels, or from some
    /// other scale detector like the Laplace scale used by Mikolajczyk
    /// & Schmid
    float scale;

    /// Integer location (unnormalized), mainly for internal use.
    int ix;
    int iy;
    //octave? plane index?

    /// Since the orientation is not necessarily unique, we may have more
    /// than one hypothesis for the orientation of an interest point.  I
    /// considered making a vector of orientations for a single point.
    /// However, it is probably better to make more than one interest
    /// point with the same (x,y,s) since the descriptor will be unique
    /// for a given orientation anyway.
    float orientation;

    /// These values are helpful for determining thresholds on the
    /// interest function value and the peak strength in the x, y and
    /// scale directions.  Types?  Maybe lT?  I'm not going to fuss at
    /// the moment because I am hoping this is a temporary thing...
    float interest;

    /// And finally the descriptor for the interest point.  SIFT points
    /// have a vector of integers, PCA-SIFT features have a vector of
    /// floats or doubles...
    vw::Vector<float> descriptor;

    int size() const { return 2; }
    double operator[] (int index) const { 
      if (index == 0) return x;
      else if (index == 1) return y;
      else vw_throw( vw::ArgumentErr() << "Interest Point: Invalid index" );
    }

    bool operator< (const InterestPoint& other) const {
      return (other.interest < interest);
    }
  };

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
  
  /// Get the size x size support region around an interest point.
  /// transform takes care of interpolation and edge extension.
  int inline get_support( ImageView<float>& support,
		   float x, float y, float scale, float ori,
		   const ImageView<float>& source, int size=SUPPORT_SIZE ) {
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

  /// Get the support region around an interest point.
  int inline get_support( ImageView<float>& support,
			  const InterestPoint& pt,
			  const ImageView<float>& source,
			  int size=SUPPORT_SIZE ) {
    return get_support(support, pt.x, pt.y, pt.scale, pt.orientation,
		       source, size);
  }

  template <class ImplT>
  struct DescriptorBase {
    /// Returns the derived implementation type.
    ImplT& impl() { return *static_cast<ImplT*>(this); }

    /// Returns the derived implementation type.
    ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

    /*
    int compute_descriptors( std::vector<InterestPoint>& points,
			     const ImageView<typename ImplT::real_type>& source ) {
      ImageView<typename ImplT::real_type> feature_support;

      impl().set_source(source);
      for (unsigned i = 0; i <points.size(); i++){
	impl().get_support(feature_support, points[i], source );
	
	// This function is specialized for different descriptors.
	impl().compute_descriptor_from_support( points[i], feature_support );
      }
      return 0;
    }

    /// Certain implementations may require support to be an array of images.
    int get_support(ImageView<typename ImplT::real_type>& support,
		    const InterestPoint& pt,
		    const ImageView<ImplT::real_type>& source) {
      vw::ip::get_support(support, pt, source); // 41x41
    }
    */
  };

  template <class DescriptorT>
  int generate_descriptors(std::vector<InterestPoint>& points,
			   const ImageView<typename DescriptorT::real_type> source,
			   DescriptorBase<DescriptorT>& desc_gen) {
    desc_gen.impl().compute_descriptors(points, source);
    return 0;
  }

  template <class T>
  class PatchDescriptor : public DescriptorBase<PatchDescriptor<T> >
  {
    ImageView<T> source;

  public:
    typedef T real_type;

    void set_source(const ImageView<T>& src) {
      source = src;
    }

    int compute_descriptor_from_support(InterestPoint& pt,
					const ImageView<T>& support) {
      pt.descriptor.set_size(SUPPORT_SIZE * SUPPORT_SIZE);
      for (int i = 0; i < SUPPORT_SIZE; i++)
	for (int j = 0; j < SUPPORT_SIZE; j++)
	  pt.descriptor(i*SUPPORT_SIZE + j) = support(i,j);
      pt.descriptor = normalize(pt.descriptor);

      return 0;
    }

    int compute_descriptors( std::vector<InterestPoint>& points,
			     const ImageView<T>& source ) {
      ImageView<T> feature_support;

      set_source(source);
      for (unsigned i = 0; i <points.size(); i++){
	get_support(feature_support, points[i], source );
	
	// This function is specialized for different descriptors.
	compute_descriptor_from_support( points[i], feature_support );
      }
      return 0;
    }

    /// Certain implementations may require support to be an array of images.
    int get_support(ImageView<T>& support,
		    const InterestPoint& pt,
		    const ImageView<T>& source) {
      vw::ip::get_support(support, pt, source); // 41x41
    }
  };

  // SaveFeatureDescriptor?

}} // namespace vw::ip

#endif //__INTERESTPOINT_DESCRIPTOR_H__
