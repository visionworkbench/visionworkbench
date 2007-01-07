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

/// \file Transform.h
/// 
/// ImageView classes for transforming the domain of an image
/// (image warping).
/// 
/// Two types of transformations are supported by these classes:
///
///
#ifndef __VW_IMAGE_TRANSFORM_H__
#define __VW_IMAGE_TRANSFORM_H__

// Vision Workbench
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/Manipulation.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

// Boost 
#include <boost/concept_check.hpp>

static const float VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE = 1;
static const float VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE = 1e10; // Ten gigapixels

namespace vw {

  // -------------------------------------------------------------------------------
  // Transform Concepts
  // -------------------------------------------------------------------------------

  /// \cond INTERNAL
  namespace transform_concepts {
    /// \endcond

    /// Transform functors are classes that provide the basic
    /// functionality defined in the concepts below.  The transform
    /// class expects mapping function to (at the very least) have
    /// defined a default constructor and a reverse() method that
    /// takes one argument: a Vector2 representing the pixel location
    /// in the output image and returns a Vector2 that is the
    /// corresponding pixel location in the input image.
    ///
    /// This concept is checked by the TransformView classes to ensure
    /// that the functor adheres to these specifications.
    template <class TransformT>
    struct TransformConcept {
      void constraints() {
        Vector2 p(0,0);
        Vector2 p_out = mapper.reverse(p);
      }
      TransformT mapper;
    };
    
    /// In addition to the requirements specified in \ref
    /// TransformConcept, the compute_transformed_bbox() functions
    /// require that the transform functor also provide the inverse of
    /// the reverse() in a method called forward().  Given a pixel
    /// location of the input image, the forward() method returns the
    /// corresponding pixel location in the output image.
    template <class TransformT>
    struct InvertibleTransformConcept {
      void constraints() {
        boost::function_requires< TransformConcept<TransformT> >();
        Vector2 p(0,0);
        Vector2 p_out = mapper.forward(p);
      }
      TransformT mapper;
    };
  }

  // -------------------------------------------------------------------------------
  // Built-in transform functors
  // -------------------------------------------------------------------------------
  template <class ImplT>
  struct TransformBase {
    /// The naive implementation.  Subclasses should override this
    /// computation with a closed form bounding box computation if one
    /// is available.
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    BBox2i compute_input_bbox(BBox2i const& output_bbox) const {
      Vector2 pt;
      BBox2i bbox;
      for( int y=output_bbox.min().y(); y<output_bbox.max().y(); ++y ) {
        for( int x=output_bbox.min().x(); x<output_bbox.max().x(); ++x ) {
          Vector2 result = impl().reverse( Vector2(x,y) );
          bbox.grow( BBox2i( (int)floor(result.x()), (int)floor(result.y()), 2, 2 ) );
        }
      }
      return bbox;
    }

    // FIXME: Should these simply be omitted?  Would that DTRT at compile time?
    inline Vector2 reverse(const Vector2 &p) const { vw_throw( NoImplErr() << "TransformBase: reverse() is not implemented for this transform function." ); return Vector2(); }
    inline Vector2 forward(const Vector2 &p) const { vw_throw( NoImplErr() << "TransformBase: forward() is not implemented for this transform function." ); return Vector2(); }
  };
  
  /// Resample Image Transform Functor
  ///
  /// Transform points for image warping by applying a (possibly
  /// non-uniform) scaling in x and y.
  class ResampleTransform : public TransformBase<ResampleTransform> {
    double m_xfactor, m_yfactor;
  public:    
    ResampleTransform(double x_scaling, double y_scaling) : 
      m_xfactor( x_scaling ) , m_yfactor( y_scaling ) {}

    inline Vector2 reverse(const Vector2 &p) const {
      return Vector2( p(0) / m_xfactor, p(1) / m_yfactor );
    }

    inline Vector2 forward(const Vector2 &p) const {
      return Vector2( p(0) * m_xfactor, p(1) * m_yfactor );
    }
  };


  /// Translate Image Transform Functor
  ///
  /// Reposition an image by applying a translation to x and y.
  class TranslateTransform : public TransformBase<TranslateTransform> {
    double m_xtrans, m_ytrans;
  public:
    TranslateTransform(double x_translation, double y_translation) : 
      m_xtrans( x_translation ) , m_ytrans( y_translation ) {}
                
    inline Vector2 reverse(const Vector2 &p) const {
      return Vector2( p(0) - m_xtrans, p(1) - m_ytrans );
    }

    inline Vector2 forward(const Vector2 &p) const {
      return Vector2( p(0) + m_xtrans, p(1) + m_ytrans );
    }
  };


  /// Homography Image Mapping Functor
  ///
  /// Transform points for image warping by applying a linear operator
  /// (a 3x3 homography).
  class HomographyTransform : public TransformBase<HomographyTransform> {
  private:
    Matrix<double> m_H_inverse;
    Matrix<double> m_H;
  public:

    HomographyTransform(Matrix<double> H) : m_H_inverse( inverse(H) ), m_H(H) {
      VW_ASSERT ( (m_H.rows() == 3) && (m_H.cols() == 3),
                  ArgumentErr() << "HomographyTransform: Invalid dimensions for homography. Matrix must be 3x3.");
    }

    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 reverse(const Vector2 &p) const {
      double w = m_H_inverse(2,0) * p(0) + m_H_inverse(2,1) * p(1) + m_H_inverse(2,2);
      return Vector2( ( m_H_inverse(0,0) * p(0) + m_H_inverse(0,1) * p(1) + m_H_inverse(0,2) ) / w,
                      ( m_H_inverse(1,0) * p(0) + m_H_inverse(1,1) * p(1) + m_H_inverse(1,2) ) / w);
    }
		
    /// This function defines the inverse transformation: from
    /// coordinates in the original image to coordinates in the
    /// transformed image.  This routine is not needed to compute the
    /// transformation, but it is used to determine the size of the
    /// output image when compute_transformed_bbox() is called.
    inline Vector2 forward(const Vector2 &p) const {
      double w = m_H(2,0) * p(0) + m_H(2,1) * p(1) + m_H(2,2);
      return Vector2( ( m_H(0,0) * p(0) + m_H(0,1) * p(1) + m_H(0,2) ) / w,
											( m_H(1,0) * p(0) + m_H(1,1) * p(1) + m_H(1,2) ) / w);
    }
  };


  /// PointLookup Image Mapping Functor
  ///
  /// Transform points in an image based on values in a lookup-table
  /// image.  The pixel location in the input image that corresponds
  /// to location (i,j) in the output image is at the position stored
  /// in lookup_image(i,j).
  class PointLookupTransform : public TransformBase<PointLookupTransform> {
    ImageView<Vector2> m_lookup_image;
  public:
    PointLookupTransform(ImageView<Vector2> &lookup_image) : m_lookup_image(lookup_image) {}
    
    inline Vector2 reverse(const Vector2 &p) const {
      VW_DEBUG_ASSERT(int(p.x()) >= 0  &&  int(p.y()) >= 0 && unsigned(p.x()) < m_lookup_image.cols() && unsigned(p.y()) < m_lookup_image.rows(),
                      LogicErr() << "Point lookup transform: exceeded lookup table dimensions.");
      return m_lookup_image(int(p.x()), int(p.y()));
    }
  };

  /// PointOffset Image Mapping Functor
  ///
  /// Transform points in an image based on offset values in a
  /// lookup-table image.  The pixel location in the input image that
  /// corresponds to location (i,j) in the output image is at the
  /// position stored in (i,j) + lookup_image(i,j).
  class PointOffsetTransform : public TransformBase<PointOffsetTransform> {
    ImageView<Vector2> m_offset_image;
  public:
    PointOffsetTransform(ImageView<Vector2> &offset_image) : m_offset_image(offset_image) {}
    
    inline Vector2 reverse(const Vector2 &p) const {
      VW_DEBUG_ASSERT(int(p.x()) >= 0  &&  int(p.y()) >= 0 && unsigned(p.x()) < m_offset_image.cols() && unsigned(p.y()) < m_offset_image.rows(),
                      LogicErr() << "Point offest transform: exceeded lookup table dimensions.");
      return p + m_offset_image(int(p.x()), int(p.y()));
    }
  };

	
  /// Radial Distortion Transform Adapter 
  /// 
  /// This is an adaptor class that allows you to write mapping
  /// functors in polar coordinates (i.e. radius, theta).  The origin
  /// of the polar coordinate system is centered at the image center,
  /// and the radius is scaled by <TT>image_width/2</TT>, so that a
  /// polar coordinate of <TT>[1.0, 0.0]</TT> places you at the center
  /// of the right edge of the image (<TT>[image_width,
  /// image_height/2]</TT> in cartesian coordinates).
  template <class ImplT, class ImageT>
  class RadialTransformAdaptor : public TransformBase< RadialTransformAdaptor<ImplT, ImageT> > {
  private:
    ImplT m_impl;
    double m_half_width;
    double m_half_height;

  public:

    RadialTransformAdaptor(ImplT const &mapping_functor, ImageViewBase<ImageT> const &image) : 
      m_impl(mapping_functor), 
      m_half_width(image.impl().cols()/2),
      m_half_height(image.impl().rows()/2) {
    }

    inline Vector2 reverse(const Vector2 &p) const {   
      boost::function_requires< transform_concepts::TransformConcept<ImplT> >();
			
      // Convert from cartesian to polar coordinates, where we scale r
      // such that r=1.0 corresponds to image_width/2, and the origin
      // is at [image_width/2, image_height/2] in the image.
      Vector2 centered_point( p(0) - m_half_width, p(1) - m_half_height );
      Vector2 polar_coordinates( norm_2(centered_point) / m_half_width, 
																 atan2(centered_point(1),centered_point(0)) );

      // Call out to the radial mapping functor that we have wrapped.
      Vector2 result = m_impl(polar_coordinates);
      
      // Convert from polar coordinates back to cartesian coordinates
      return Vector2(m_half_width * result(0) * cos(result(1)) + m_half_width, 
										 m_half_width * result(0) * sin(result(1)) + m_half_height);
    }

    inline Vector2 forward(const Vector2 &p) const {
      boost::function_requires< transform_concepts::InvertibleTransformConcept<ImplT> >();
     
      Vector2 centered_point( p(0) - m_half_width, p(1) - m_half_height );
      Vector2 polar_coordinates(norm_2(centered_point) / m_half_width, 
																atan2(centered_point(1),centered_point(0)) );
      Vector2 result = m_impl.reverse(polar_coordinates);
      return Vector2(m_half_width * result(0) * cos(result(1)) + m_half_width, 
										 result(0) * sin(result(1)) + m_half_height);
    }
  };
	

  /// Transform Composition Transform Adapter 
  /// 
  /// This is a wrapper class that allows you to composte two transform 
  /// functors.  The arguments to the constructor are in the usual 
  /// function composition order.  That is, CompositionTransform(tx1,tx2)
  /// yields a functor whose forward mapping is tx1.forward(tx2.forward(v)), 
  /// and thus whose reverse mapping is tx2.reverse(tx1.reverse(v)).  
  /// The usual way to construct a CompositionTransform is with the 
  /// compose() free function, below.
  template <class Tx1T, class Tx2T>
  class CompositionTransform : public TransformBase<CompositionTransform<Tx1T,Tx2T> >
  {
    Tx1T tx1;
    Tx2T tx2;
  public:
    CompositionTransform( Tx1T const& tx1, Tx2T const& tx2 ) : tx1(tx1), tx2(tx2) {}
    
    inline Vector2 forward( Vector2 const& p ) const {
      return tx1.forward( tx2.forward( p ) );
    }
    
    inline Vector2 reverse( Vector2 const& p ) const {
      return tx2.reverse( tx1.reverse( p ) );
    }
  };

  /// Composes two transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T>
  CompositionTransform<Tx1T,Tx2T>
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2 ) {
    typedef CompositionTransform<Tx1T,Tx2T> result_type;
    return result_type( tx1.impl(), tx2.impl() );
  }

  /// Composes three transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T, class Tx3T>
  CompositionTransform<Tx1T,CompositionTransform<Tx2T,Tx3T> > 
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2,
                  TransformBase<Tx3T> const& tx3 ) {
    typedef CompositionTransform<Tx1T,CompositionTransform<Tx2T,Tx3T> > result_type;
    return result_type( tx1.impl(), compose( tx2, tx3 ) );
  }

  /// Composes three transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T, class Tx3T, class Tx4T>
  CompositionTransform<Tx1T,CompositionTransform<Tx2T,CompositionTransform<Tx3T,Tx4T> > >
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2,
                  TransformBase<Tx3T> const& tx3,
                  TransformBase<Tx3T> const& tx4 ) {
    typedef CompositionTransform<Tx1T,CompositionTransform<Tx2T,CompositionTransform<Tx3T,Tx4T> > > result_type;
    return result_type( tx1.impl(), compose( tx2, tx3, tx4 ) );
  }


  // ------------------------
  // class TransformView
  // ------------------------

  /// An image view for transforming an image with an arbitrary mapping functor.  
  template <class ImageT, class TransformT>
  class TransformView : public ImageViewBase<TransformView<ImageT,TransformT> >
  {

    ImageT m_image;
    TransformT m_mapper;
    unsigned m_width, m_height;
    
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<TransformView> pixel_accessor;

    /// \cond INTERNAL
    BOOST_CLASS_REQUIRE(TransformT, transform_concepts, TransformConcept);
    /// \endcond
    
    /// The default constructor creates a tranformed image with the
    /// same dimensions as the original.
    TransformView( ImageT const& view, TransformT const& mapper ) :
      m_image(view), m_mapper(mapper), m_width(view.cols()), m_height(view.rows()) {}

    /// This constructor allows you to specify the size of the
    /// transformed image.
    TransformView( ImageT const& view, TransformT const& mapper, unsigned width, unsigned height ) : 
      m_image(view), m_mapper(mapper), m_width(width), m_height(height) {}
    
    inline unsigned cols() const { return m_width; }
    inline unsigned rows() const { return m_height; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline result_type operator()( float i, float j, int p=0 ) const {
      Vector2 pt_backprojected = m_mapper.reverse(Vector2(i,j));
      return m_image(pt_backprojected[0], pt_backprojected[1], p);
    }

    /// \cond INTERNAL
    typedef TransformView<typename ImageT::prerasterize_type, TransformT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox) const { 
      BBox2i transformed_bbox = m_mapper.compute_input_bbox(bbox);
      return prerasterize_type( m_image.prerasterize(transformed_bbox), m_mapper, m_width, m_height );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { 
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type Traits
  template <class ImplT, class TransformT>
  struct IsFloatingPointIndexable<TransformView<ImplT, TransformT> > : public true_type {};
  /// \endcond



  /// Compute the bounding box in the transformed image space that
  /// contains all of the transformed pixels.  The bounding box is
  /// computed by forward transforming all of the pixel coordinates
  /// from the original image.  If you know your transformation to be
  /// convex, you will probably want to use
  /// compute_transformed_bbox_fast(), below.  If the bounding box
  /// exceeds preset limits, a warning will be printed out.
  template <class ViewT, class TransformT>
  inline BBox2f compute_transformed_bbox(ImageViewBase<ViewT> const& image, 
                                         TransformT const& transform_func,
                                         float min_image_size = VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE,
                                         float max_image_size = VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE) {
    Vector2 pt;
    BBox2f bbox;
    for (pt[0] = 0; pt[0] < image.impl().cols(); (pt[0])++)
      for (pt[1] = 0; pt[1] < image.impl().rows(); (pt[1])++) 
        bbox.grow(transform_func.forward(pt));
    
    // If the image bounding box is too large or too small, print a
    // warning message.
    if ( (bbox.width() * bbox.height()) < min_image_size ) {
      std::cout << "Warning: The transformed image exceeds the minimum (" << min_image_size << ") \n"
                << "         recommended image dimension in compute_transformed_bbox().\n";
    }
    if ( (bbox.width() * bbox.height()) > max_image_size ) {
      std::cout << "Warning: The transformed image exceeds the maximum (" << max_image_size << ") \n"
                << "         recommended image dimension in compute_transformed_bbox().\n";
    }
    return bbox;
  }

  /// Compute full extent of the transformed image by performing the
  /// inverse transformation on the pixel coordinates of the input
  /// image.  This version of the function computes the bounds by
  /// computing only the transformed pixel locations for the perimiter
  /// pixels of the input image.  This is much faster than calculating
  /// the transformed position for all pixels in the input image, and
  /// for most transformations points on the interior of the input
  /// image will end up on the interior of the output image.
  template <class ViewT, class TransformT>
  inline BBox2f compute_transformed_bbox_fast(ImageViewBase<ViewT> const& image, 
                                              TransformT const& transform_func,
                                              unsigned min_image_size = VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE,
                                              unsigned max_image_size = VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE) {
    Vector2 pt;
    BBox2f bbox;
    
    // Top edge
    for (pt[0] = 0; pt[0] < image.impl().cols(); (pt[0])++) {
      pt[1] = 0;
      bbox.grow(transform_func.forward(pt));
    }

    // Bottom edge
    for (pt[0] = 0; pt[0] < image.impl().cols(); (pt[0])++) {
      pt[1] = image.impl().rows() - 1;
      bbox.grow(transform_func.forward(pt));
    }

    // Left edge
    for (pt[1] = 0; pt[1] < image.impl().rows(); (pt[1])++) {
      pt[0] = 0;
      bbox.grow(transform_func.forward(pt));
    }

    // Right edge
    for (pt[1] = 0; pt[1] < image.impl().rows(); (pt[1])++) {
      pt[0] = image.impl().cols() - 1;
      bbox.grow(transform_func.forward(pt));
    }
    
    // If the image bounding box is too large or too small, fall
    // back and set the size of the output image to the size of the
    // input image and print a warning message.
    if ( (bbox.width() * bbox.height()) < min_image_size ) {
      std::cout << "Warning: The transformed image exceeds the minimum (" << min_image_size << ") \n"
                << "         recommended image dimension in compute_transformed_bbox_fast().\n";
    }
    if ( (bbox.width() * bbox.height()) > max_image_size ) {
      std::cout << "Warning: The transformed image exceeds the maximum (" << max_image_size << ") \n"
                << "         recommended image dimension in compute_transformed_bbox_fast().\n";
    }
    return bbox;
  }

  // -------------------------------------------------------------------------------
  // Functional API
  // -------------------------------------------------------------------------------

  /// Apply a Transformation using an arbitary mapping functor.
  ///
  /// Returns a transformed image view.  The user can specify the type
  /// of interpolation and edge extension to be done by supplying the
  /// appropriate functors in the last two arguments.  For example:
  /// 
  /// <TT>transform(input, transform_functor, BicubicInterpolation(), ZeroEdgeExtension());</TT>
  ///
  /// See Interpolation.h and EdegExtension.h for a list of built-in
  /// interpolation and edge extension modes.  It is also possible to
  /// write your own.  Again, see the above files for more information.
  template <class ImageT, class TransformT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT> >::type
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func, 
                    EdgeT const& edge_func, InterpT const& interp_func) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
      (interpolate(v, interp_func, edge_func), transform_func);
  }
  
  /// Convenience function: transform with a default interpolation
  /// scheme of bilinear interpolation.
  template <class ImageT, class TransformT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func, 
                    EdgeT const& edge_func) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT> 
      (interpolate(v, BilinearInterpolation(), edge_func), transform_func);
  }
  
  /// Convenience function: transform with a default scheme of
  /// bilinear interpolation and zero edge extension.
  template <class ImageT, class TransformT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT> 
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func);
  }





  /// This variant of transform allows the user to specify the
  /// dimensions of the transformed image.  The upper left hand point
  /// (0,0) stays fixed.  For a more flexible method of cropping to an
  /// arbitrary bounding box, use one of the transform methods defined
  /// below.
  template <class ImageT, class TransformT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func, 
                    int width, int height, EdgeT const& edge_func, InterpT const& interp_func) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
      (interpolate(v, interp_func, edge_func), transform_func, width, height);
  }
  
  /// Convenience function: transform with a default interpolation
  /// scheme of bilinear interpolation. The user can specify the
  /// dimensions of the output image.
  template <class ImageT, class TransformT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func, 
                    int width, int height, EdgeT const& edge_func) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT> 
      (interpolate(v, BilinearInterpolation(), edge_func), transform_func, width, height);
  }

  /// Convenience function: transform with a default scheme of
  /// bilinear interpolation and zero edge extension.  The user can
  /// specify the dimensions of the output image.
  template <class ImageT, class TransformT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v, TransformT const& transform_func,
                    int width, int height) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT> 
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func, width, height);
  }
  





  /// Transform an image and select the bounding box in the
  /// transformed space from which to take the new image.  The
  /// compute_transformed_bbox() method can be used to compute the
  /// bounding box that will fit all of the transformed pixels.
  template <class ImageT, class TransformT, class BBoxRealT, class EdgeT, class InterpT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v, 
                    TransformT const& transform_func, 
                    BBox<BBoxRealT,2> const& bbox, 
                    EdgeT const& edge_func,
                    InterpT const& interp_func) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
                (interpolate(v, interp_func, edge_func), transform_func), bbox);
  }

  /// Convenience function: free transform with a default scheme of
  /// bilinear interpolation.  The user supplies a bounding box in the
  /// transformed space from that determines what pixels will be
  /// rasterized.
  template <class ImageT, class TransformT, class BBoxRealT, class EdgeT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v, 
                    TransformT const& transform_func, 
                    BBox<BBoxRealT,2> const& bbox, 
                    EdgeT const& edge_func) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
                (interpolate(v, BilinearInterpolation(), edge_func), transform_func), bbox);
  }

  /// Convenience function: transform with a default scheme of
  /// bilinear interpolation and zero edge extension. The user
  /// supplies a bounding box in the transformed space from that
  /// determines what pixels will be rasterized.
  template <class ImageT, class TransformT, class BBoxRealT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v, 
                    TransformT const& transform_func, 
                    BBox<BBoxRealT,2> const& bbox) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
                (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func), bbox);
  }



  // -------------------------------------------------------------------------------
  // Resample
  // -------------------------------------------------------------------------------

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v, double x_scale_factor, double y_scale_factor, 
                   EdgeT const& edge_func, InterpT const& interp_func) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor), 
                     int(round(v.impl().cols()*x_scale_factor)), int(round(v.impl().rows()*y_scale_factor)),
                     edge_func, interp_func);
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, ResampleTransform> 
  inline resample( ImageViewBase<ImageT> const& v, double x_scale_factor, double y_scale_factor, 
                   EdgeT const& edge_func) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor), 
                     int(round(v.impl().cols()*x_scale_factor)), int(round(v.impl().rows()*y_scale_factor)),
                     edge_func, vw::BilinearInterpolation());
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform> 
  inline resample( ImageViewBase<ImageT> const& v, 
                   double x_scale_factor, double y_scale_factor) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor), 
                     int(round(v.impl().cols()*x_scale_factor)), int(round(v.impl().rows()*y_scale_factor)), 
                     vw::ConstantEdgeExtension(), vw::BilinearInterpolation());
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y and the dimensions of the output image.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v, double scale_factor, 
                   int output_width, int output_height,
                   EdgeT const& edge_func, InterpT const& interp_func) {
    return transform(v, ResampleTransform(scale_factor, scale_factor), 
                     output_width, output_height, 
                     edge_func, interp_func);
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y and the dimensions of the output image.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v, double scale_factor, 
                   EdgeT const& edge_func, InterpT const& interp_func) {
    return transform(v, ResampleTransform(scale_factor, scale_factor), 
                     int(round(v.impl().cols()/scale_factor)), int(round(v.impl().rows()*scale_factor)),
                     edge_func, interp_func);
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y and the dimensions of the output image.
  template <class ImageT, class EdgeT>
  typename boost::disable_if<IsScalar<EdgeT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, ResampleTransform> >::type 
  inline resample( ImageViewBase<ImageT> const& v, double scale_factor, EdgeT const& edge_func) {
    return transform(v, ResampleTransform(scale_factor, scale_factor),
                     int(round(v.impl().cols()*scale_factor)), int(round(v.impl().rows()*scale_factor)),
                     edge_func, vw::BilinearInterpolation());
  }

  /// Resample the image.  The user specifies the scaling factor in x
  /// and y and the dimensions of the output image.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform> 
  inline resample( ImageViewBase<ImageT> const& v, double scale_factor) {
    return transform(v, ResampleTransform(scale_factor, scale_factor), 
                     int(round(v.impl().cols()*scale_factor)), int(round(v.impl().rows()*scale_factor)),
                     vw::ConstantEdgeExtension(), vw::BilinearInterpolation());
  }


  // -------------------------------------------------------------------------------
  // Translate
  // -------------------------------------------------------------------------------

  /// Translate the image.  The user specifies the offset in x
  /// and y.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TranslateTransform> >::type
  inline translate( ImageViewBase<ImageT> const& v, double x_offset, double y_offset,
                   EdgeT const& edge_func, InterpT const& interp_func) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     edge_func, interp_func);
  }

  /// Translate the image.  The user specifies the offset in x
  /// and y.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TranslateTransform> 
  inline translate( ImageViewBase<ImageT> const& v, double x_offset, double y_offset, 
                   EdgeT const& edge_func) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     edge_func, BilinearInterpolation());
  }

  /// Translate the image.  The user specifies the offset in x
  /// and y.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TranslateTransform> 
  inline translate( ImageViewBase<ImageT> const& v, 
                   double x_offset, double y_offset) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     ZeroEdgeExtension(), BilinearInterpolation());
  }

} // namespace vw

#endif // __VW_IMAGE_TRANSFORM_H__ 
