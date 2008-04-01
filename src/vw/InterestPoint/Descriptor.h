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
#ifndef __VW_INTERESTPOINT_DESCRIPTOR_H__
#define __VW_INTERESTPOINT_DESCRIPTOR_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/InterestPoint/InterestData.h>

#include <vw/InterestPoint/MatrixIO.h>
#include <vw/InterestPoint/VectorIO.h>

namespace vw { 
namespace ip {

  template <class ImplT>
  class DescriptorGeneratorBase {

    static const int DEFAULT_SUPPORT_SIZE = 41;

  public:

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond


    // Given an image and a list of interest points, set the
    // descriptor field of the interest points using the
    // compute_descriptor() method provided by the subclass.
    template <class ViewT>
    void operator() ( ImageViewBase<ViewT> const& image, InterestPointList& points ) {
      //vw_out(InfoMessage) << "\nComputing support and descriptor for " << points.size() << " interest points\n";
      for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {

        // First we compute the support region based on the interest point 
        ImageView<PixelGray<float> > support = get_support(*i, pixel_cast<PixelGray<float> >(channel_cast_rescale<float>(image.impl())));

        // Next, we pass the support region and the interest point to
        // the descriptor generator ( compute_descriptor() ) supplied
        // by the subclass.
        i->descriptor = impl().compute_descriptor(support);
      }
    }

    /// Get the size x size support region around an interest point,
    /// rescaled by the scale factor and rotated by the specified
    /// angle.
    template <class ViewT>
    inline ImageView<typename ViewT::pixel_type> get_support( float x, float y, float scale, float ori,
                                                              ImageViewBase<ViewT> const& source, int size=DEFAULT_SUPPORT_SIZE ) {
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
    get_support( InterestPoint const& pt, ImageViewBase<ViewT> const& source, int size=DEFAULT_SUPPORT_SIZE ) {
      return get_support(pt.x, pt.y, pt.scale, pt.orientation, source.impl(), size);
    }
    
  };



  /// A basic example descriptor class. The descriptor for an interest
  /// point is simply the pixel values in the support region around
  /// the point. It is normalized to provide some tolerance to changes
  /// in illumination.
  struct PatchDescriptorGenerator : public DescriptorGeneratorBase<PatchDescriptorGenerator> {

    template <class ViewT>
    Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) const {
      Vector<float> result;
      result.set_size(support.impl().cols() * support.impl().rows());

      for (int j = 0; j < support.impl().rows(); j++)
        for (int i = 0; i < support.impl().cols(); i++) {
          PixelGray<float> pix(support.impl()(i,j));
          result(j*support.impl().cols() + i) = pix.v();
        }

      return normalize(result);
    }

  };

  // An implementation of PCA-SIFT
  struct PCASIFTDescriptorGenerator : public DescriptorGeneratorBase<PCASIFTDescriptorGenerator> {

    std::string basis_filename, avg_filename;
    
    PCASIFTDescriptorGenerator(const std::string& pcabasis_filename,
			       const std::string& pcaavg_filename) 
      : basis_filename(pcabasis_filename), avg_filename(pcaavg_filename) { }


    template <class ViewT>
    Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) const {
      Matrix<float> pca_basis;
      Vector<float> pca_avg;
      
      // TODO: Move read operation out to the constructor
      read_matrix(pca_basis, basis_filename);
      read_vector(pca_avg, avg_filename);

      Vector<float> result;
      result.set_size(pca_basis.cols());

      // compute normalization constant (sum squares)
      double norm_const = 0;
      for (int j = 0; j < support.impl().rows(); j++) {
	for (int i = 0; i < support.impl().cols(); i++) {
	  norm_const += support.impl()(i,j) * support.impl()(i,j);
	}
      }
      norm_const = sqrt(norm_const);

      double norm_pixel = 0;
      // project image patch onto PCA basis to get descriptor
      unsigned int index = 0;
      for (int j = 0; j < support.impl().rows(); j++) {
        for (int i = 0; i < support.impl().cols(); i++) {
	  norm_pixel = support.impl()(i,j).v()/norm_const - pca_avg(index);
	  
	  for (int k = 0; k < pca_basis.cols(); k++) {
	    result[k] += norm_pixel * pca_basis(index,k);
	  }
	  ++index;
	}
      }
      return result;
    }
  };
  
}} // namespace vw::ip


#endif //__VW_INTERESTPOINT_DESCRIPTOR_H__
