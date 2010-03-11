// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Descriptor.h
///
/// Basic classes and functions for generating interest point descriptors.
///
#ifndef __VW_INTERESTPOINT_DESCRIPTOR_H__
#define __VW_INTERESTPOINT_DESCRIPTOR_H__

#include <vw/Core/Debugging.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/MatrixIO.h>
#include <vw/InterestPoint/VectorIO.h>
#include <vw/InterestPoint/IntegralImage.h>

namespace vw {
namespace ip {

  template <class ImplT>
  class DescriptorGeneratorBase {

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

      // Timing
      Timer total("\tTotal elapsed time", DebugMessage, "interest_point");

      for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {

        // First we compute the support region based on the interest point
        ImageView<PixelGray<float> > support = get_support(*i, pixel_cast<PixelGray<float> >(channel_cast_rescale<float>(image.impl())),
                                                           impl().support_size() );

        // Next, we pass the support region and the interest point to
        // the descriptor generator ( compute_descriptor() ) supplied
        // by the subclass.
        i->descriptor = impl().compute_descriptor(support);
      }

    }

    // Default suport size ( i.e. descriptor window)
    int support_size( void ) {
      return 41;
    }

    /// Get the size x size support region around an interest point,
    /// rescaled by the scale factor and rotated by the specified
    /// angle.
    template <class ViewT>
    inline ImageView<typename ViewT::pixel_type> get_support( float x, float y, float scale, float ori,
                                                              ImageViewBase<ViewT> const& source, int size) {

      float half_size = ((float)(size - 1)) / 2.0f;
      float scaling = 1.0f / scale;

      return transform(source.impl(),
                       compose(TranslateTransform(half_size, half_size),
                               compose(ResampleTransform(scaling, scaling),
                                       RotateTransform(-ori),
                                       TranslateTransform(-x, -y))),
                       size, size );
    }

    /// Get the support region around an interest point, scaled and
    /// rotated appropriately.
    template <class ViewT>
    inline ImageView<typename ViewT::pixel_type>
    get_support( InterestPoint const& pt, ImageViewBase<ViewT> const& source, int size ) {
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
    Matrix<float> pca_basis;
    Vector<float> pca_avg;

    PCASIFTDescriptorGenerator(const std::string& pcabasis_filename,
                               const std::string& pcaavg_filename)
      : basis_filename(pcabasis_filename), avg_filename(pcaavg_filename) {

      // Read the PCA basis matrix and average vector
      read_matrix(pca_basis, basis_filename);
      read_vector(pca_avg, avg_filename);
    }


    template <class ViewT>
    Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) const {

      Vector<float> result(pca_basis.cols());

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

          for (unsigned k = 0; k < pca_basis.cols(); k++) {
            result[k] += norm_pixel * pca_basis(index,k);
          }
          ++index;
        }
      }
      return result;
    }
  };

  // A Simple Scaled Gradient descriptor that reduces the number of elements
  // used in the descriptor and is hopefully more robust against
  // illumination changes.
  struct SGradDescriptorGenerator : public DescriptorGeneratorBase<SGradDescriptorGenerator> {

    static const uint box_strt[5];
    static const uint box_size[5];
    static const uint box_half[5];

    template <class ViewT>
    Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) const {
      Vector<float> result(180);

      typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;
      ImageView<channel_type> iimage = IntegralImage(support);

      uint write_idx = 0;

      // Iterate through scales
      for ( uint s = 0; s < 5; s++ ) {
        float inv_bh2 = 1 / float(box_half[s]*box_half[s]);

        // Iterate though quadrants
        for ( uint i = 0; i < 3; i++ ) {
          for ( uint j = 0; j < 3; j++ ) {
            Vector2i top_left( box_strt[s]+i*box_size[s],
                               box_strt[s]+j*box_size[s] );

            float minor_quad[4];

            // 1.) Top Left in local
            minor_quad[0] = IntegralBlock( iimage,
                                           top_left,
                                           top_left+Vector2i(box_half[s],
                                                             box_half[s]) );
            // 2.) Top Right in local
            minor_quad[1] = IntegralBlock( iimage,
                                           top_left+Vector2i(box_half[s],0),
                                           top_left+Vector2i(2*box_half[s],
                                                             box_half[s]) );
            // 3.) Bot Left in local
            minor_quad[2] = IntegralBlock( iimage,
                                           top_left+Vector2i(0,box_half[s]),
                                           top_left+Vector2i(box_half[s],
                                                             2*box_half[s]) );
            // 4.) Bot Right in local
            minor_quad[3] = IntegralBlock( iimage,
                                           top_left+Vector2i(box_half[s],
                                                             box_half[s]),
                                           top_left+Vector2i(2*box_half[s],
                                                             2*box_half[s]) );

            // 5.) Pulling out gradients
            result[write_idx] = (minor_quad[0] - minor_quad[2])*inv_bh2; write_idx++;
            result[write_idx] = (minor_quad[1] - minor_quad[3])*inv_bh2; write_idx++;
            result[write_idx] = (minor_quad[0] - minor_quad[1])*inv_bh2; write_idx++;
            result[write_idx] = (minor_quad[2] - minor_quad[3])*inv_bh2; write_idx++;
          } // end j
        } // end i
      } // end s
      return normalize(result);
    }

    int support_size( void ) { return 42; }
  };

}} // namespace vw::ip


#endif //__VW_INTERESTPOINT_DESCRIPTOR_H__
