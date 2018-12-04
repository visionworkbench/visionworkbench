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


// \file Transform.h
//
// ImageView classes for transforming the domain of an image
// (image warping).
//
//
#ifndef __VW_IMAGE_TRANSFORM_H__
#define __VW_IMAGE_TRANSFORM_H__

#include <vw/Math/Transform.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>

static const double VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE = 1;
static const double VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE = 1e10; // Ten gigapixels

namespace vw {

  // PointLookup image transform functor
  //
  // Transform points in an image based on values in a lookup-table
  // image.  The pixel location in the input image that corresponds
  // to location (i,j) in the output image is at the position stored
  // in lookup_image(i,j).
  class PointLookupTransform : public TransformBase<PointLookupTransform> {
    ImageViewRef<Vector2> m_lookup_image;
  public:
    template <class OffsetImageT>
    PointLookupTransform(ImageViewBase<OffsetImageT> const& lookup_image)
      : m_lookup_image(lookup_image.impl()) {}

    inline Vector2 reverse(const Vector2 &p) const {
      int32 x = (int32) p.x(), y = (int32) p.y();
      VW_DEBUG_ASSERT(x>=0 && y>=0 && x<m_lookup_image.cols() && y<m_lookup_image.rows(),
                      LogicErr() << "Point lookup transform: exceeded lookup table dimensions.");
      return m_lookup_image(x,y);
    }
  };


  // PointOffset image transform functor
  //
  // Transform points in an image based on offset values in a
  // lookup-table image.  The pixel location in the input image that
  // corresponds to location (i,j) in the output image is at the
  // position stored in (i,j) + lookup_image(i,j).
  class PointOffsetTransform : public TransformBase<PointOffsetTransform> {
    ImageViewRef<Vector2> m_offset_image;
  public:
    template <class OffsetImageT>
    PointOffsetTransform(ImageViewBase<OffsetImageT> const& offset_image)
      : m_offset_image(offset_image.impl()) {}

    inline Vector2 reverse(const Vector2 &p) const {
      int32 x = (int32) p.x(), y = (int32) p.y();
      VW_DEBUG_ASSERT(x>=0 && y>=0 && x<m_offset_image.cols() && y<m_offset_image.rows(),
                      LogicErr() << "Point offest transform: exceeded lookup table dimensions.");
      return p + m_offset_image(x,y);
    }
  };


  // RadialTransformAdaptor image transform functor adaptor
  //
  // This is an adaptor class that allows you to write mapping
  // functors in polar coordinates (i.e. radius, theta).  The origin
  // of the polar coordinate system is centered at the image center,
  // and the radius is scaled by <TT>image_width/2</TT>, so that a
  // polar coordinate of <TT>[1.0, 0.0]</TT> places you at the center
  // of the right edge of the image (<TT>[image_width,
  // image_height/2]</TT> in cartesian coordinates).
  template <class ImplT, class ImageT>
  class RadialTransformAdaptor : public TransformBase<RadialTransformAdaptor<ImplT, ImageT> > {
  private:
    ImplT  m_impl;
    double m_half_width;
    double m_half_height;

  public:

    RadialTransformAdaptor(ImplT const &mapping_functor, ImageViewBase<ImageT> const &image) :
      m_impl(mapping_functor),
      m_half_width (image.impl().cols()/2),
      m_half_height(image.impl().rows()/2) {
    }

    inline Vector2 reverse(const Vector2 &p) const {
      // Convert from cartesian to polar coordinates, where we scale r
      // such that r=1.0 corresponds to image_width/2, and the origin
      // is at [image_width/2, image_height/2] in the image.
      Vector2 centered_point   ( p(0) - m_half_width, p(1) - m_half_height );
      Vector2 polar_coordinates( norm_2(centered_point) / m_half_width,
                                 atan2(centered_point(1),centered_point(0)) );

      // Call out to the radial mapping functor that we have wrapped.
      Vector2 result = m_impl(polar_coordinates);

      // Convert from polar coordinates back to cartesian coordinates
      return Vector2(m_half_width * result(0) * cos(result(1)) + m_half_width,
                     m_half_width * result(0) * sin(result(1)) + m_half_height);
    }

    inline Vector2 forward(const Vector2 &p) const {
      Vector2 centered_point   ( p(0) - m_half_width, p(1) - m_half_height );
      Vector2 polar_coordinates( norm_2(centered_point) / m_half_width,
                                 atan2(centered_point(1),centered_point(0)) );
      Vector2 result = m_impl.reverse(polar_coordinates);
      return Vector2(m_half_width * result(0) * cos(result(1)) + m_half_width,
                     result(0) * sin(result(1)) + m_half_height);
    }

    virtual FunctionType forward_type() const { return m_impl.forward_type(); }
    virtual FunctionType reverse_type() const { return m_impl.reverse_type(); }
  };

  // ApproximateTransform image transform functor template.
  //
  // Mimics the behavior of a given transform functor, but attempts to
  // build a lookup table to linearly interpolate approimate results
  // to the reverse() function for arguments within the given bounding
  // box, to within the original transform functor's tolerance.
  template <class TransformT>
  class ApproximateTransform : public TransformT {
    BBox2i m_bbox;
    ImageView<Vector2> m_table;
  public:
    ApproximateTransform( TransformT const& transform, BBox2i const& bbox )
      : TransformT( transform ), m_bbox( bbox ), m_table(2,2)
    {
      // Initialize with a simple 2x2 lookup table
      int32 n=2;
      m_table(0,0) = TransformT::reverse(bbox.min());
      m_table(1,0) = TransformT::reverse(Vector2(bbox.max().x(),bbox.min().y()));
      m_table(0,1) = TransformT::reverse(Vector2(bbox.min().x(),bbox.max().y()));
      m_table(1,1) = TransformT::reverse(bbox.max());

      // Double the grid density until the worst (squared) approximation error
      // is less than the allowed (squared) tolerance.
      double max_sqr_err = 0;
      double tol_sqr = TransformT::tolerance() * TransformT::tolerance();
      Vector2 origin = bbox.min(), diag = bbox.size();
      do {
        n = 2*n-1;
        // Fall back for unapproximatably crazy transform functions.
        if( n>=bbox.width()|| n>=bbox.height() ) {
          m_table.reset();
          return;
        }
        ImageView<Vector2> prev = m_table;
        m_table.set_size(n,n);
        max_sqr_err = 0;
        for( int y=0; y<n; ++y ) {
          for( int x=0; x<n; ++x ) {
            if( (y%2)==0 && (x%2==0) ) {
              m_table(x,y) = prev(x/2,y/2);
            }
            else {
              Vector2 pos = Vector2(x,y)/(n-1);
              m_table(x,y) = TransformT::reverse(origin+elem_prod(pos,diag));
              Vector2 interp;
              if( (y%2)==0 ) interp = (prev(x/2,y/2) + prev(x/2+1,y/2)) / 2.0;
              else if( (x%2)==0 ) interp = (prev(x/2,y/2) + prev(x/2,y/2+1)) / 2.0;
              else interp = (prev(x/2,y/2) + prev(x/2,y/2+1) + prev(x/2+1,y/2) + prev(x/2+1,y/2+1)) / 4.0;
              double sqr_err = norm_2_sqr( m_table(x,y) - interp );
              if( sqr_err > max_sqr_err ) max_sqr_err = sqr_err;
            }
          }
        }
      } while( max_sqr_err > tol_sqr );
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      // Fall back if the function was not approximatable.
      if( ! m_table.is_valid_image() )
        return TransformT::reverse( p );

      // We re-implement bilinear interpolation by hand here because for
      // some reason the BilinearInterpolation object is exceptionally
      // slow for Vector data still.
      int    n  = m_table.cols() - 1;
      double px = n * (p.x() - m_bbox.min().x()) / (m_bbox.max().x() - m_bbox.min().x());
      double py = n * (p.y() - m_bbox.min().y()) / (m_bbox.max().y() - m_bbox.min().y());
      int32  ix = math::impl::_floor(px);
      if( ix < 0  ) ix = 0;
      if( ix >= n ) ix = n-1;
      int32 iy = math::impl::_floor(py);
      if( iy < 0  ) iy = 0;
      if( iy >= n ) iy = n-1;
      double normx = px-ix, normy = py-iy;

      Vector2 const& m00 = m_table(ix,  iy);
      Vector2 const& m10 = m_table(ix+1,iy);
      Vector2 const& m01 = m_table(ix,  iy+1);
      Vector2 const& m11 = m_table(ix+1,iy+1);

      return Vector2( (m00.x()*(1-normy)+m01.x()*normy)*(1-normx) +
                      (m10.x()*(1-normy)+m11.x()*normy)*normx,
                      (m00.y()*(1-normy)+m01.y()*normy)*(1-normx) +
                      (m10.y()*(1-normy)+m11.y()*normy)*normx );
    }

    // Never re-approximate the approximation.
    virtual double tolerance() const { return 0; }

  };

  // ------------------------
  // compute_transformed_bbox functions
  // - These could maybe go to /Math/Transforms.h
  // ------------------------

  // Compute the bounding box in the transformed image space that
  // contains all of the transformed pixels.  The bounding box is
  // computed by forward transforming all of the pixel coordinates
  // from the original image.  If you know your transformation to be
  // convex, you will probably want to use
  // compute_transformed_bbox_fast(), below.  If the bounding box
  // exceeds preset limits, a warning will be printed out.
  template <class TransformT>
  inline BBox2f compute_transformed_bbox(Vector2i   const& image_size,
                                         TransformT const& transform_func,
                                         double min_image_size = VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE,
                                         double max_image_size = VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE) {
    Vector2 pt;
    BBox2f  bbox;
    for (pt[0] = 0; pt[0] < image_size[0]; (pt[0])++)
      for (pt[1] = 0; pt[1] < image_size[1]; (pt[1])++)
        bbox.grow(transform_func.forward(pt));

    // If the image bounding box is too large or too small, print a
    // warning message.
    if ( (bbox.width() * bbox.height()) < min_image_size ) {
      VW_OUT(WarningMessage, "image") << "Warning: The transformed image exceeds the minimum (" << min_image_size << ") \n"
                                      << "         recommended image dimension in compute_transformed_bbox().\n";
    }
    if ( (bbox.width() * bbox.height()) > max_image_size ) {
      VW_OUT(WarningMessage, "image") << "Warning: The transformed image exceeds the maximum (" << max_image_size << ") \n"
                                      << "         recommended image dimension in compute_transformed_bbox().\n";
    }
    return bbox;
  }


  /// Compute full extent of the transformed image by performing the
  /// inverse transformation on the pixel coordinates of the input
  /// image.  This version of the function computes the bounds by
  /// computing only the transformed pixel locations for the perimiter
  /// pixels of the ROI in the input image.  This is much faster than calculating
  /// the transformed position for all pixels in the input image, and
  /// for most transformations points on the interior of the input
  /// image will end up on the interior of the output image.
  template <class TransformT>
  inline BBox2f compute_transformed_bbox_fast(BBox2i     const& image_roi,
                                              TransformT const& transform_func,
                                              double min_image_size = VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE,
                                              double max_image_size = VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE) {
    Vector2 pt;
    BBox2f  bbox;

    // Top edge
    for (pt[0] = image_roi.min()[0]; pt[0] < image_roi.max()[0]; (pt[0])++) {
      pt[1] = image_roi.min()[1];
      bbox.grow(transform_func.forward(pt));
    }

    // Bottom edge
    for (pt[0] = image_roi.min()[0]; pt[0] < image_roi.max()[0]; (pt[0])++) {
      pt[1] = image_roi.max()[1] - 1;
      bbox.grow(transform_func.forward(pt));
    }

    // Left edge
    for (pt[1] = image_roi.min()[1]; pt[1] < image_roi.max()[1]; (pt[1])++) {
      pt[0] = image_roi.min()[0];
      bbox.grow(transform_func.forward(pt));
    }

    // Right edge
    for (pt[1] = image_roi.min()[1]; pt[1] < image_roi.max()[1]; (pt[1])++) {
      pt[0] = image_roi.max()[0] - 1;
      bbox.grow(transform_func.forward(pt));
    }

    // If the image bounding box is too large or too small, fall
    // back and set the size of the output image to the size of the
    // input image and print a warning message.
    if ( (bbox.width() * bbox.height()) < min_image_size ) {
      VW_OUT(WarningMessage, "image") << "Warning: The transformed image exceeds the minimum (" << min_image_size << ") \n"
                                      << "         recommended image dimension in compute_transformed_bbox_fast().\n";
    }
    if ( (bbox.width() * bbox.height()) > max_image_size ) {
      VW_OUT(WarningMessage, "image") << "Warning: The transformed image exceeds the maximum (" << max_image_size << ") \n"
                                      << "         recommended image dimension in compute_transformed_bbox_fast().\n";
    }
    return bbox;
  }

  /// Overload that starts the ROI at 0,0
  template <class TransformT>
  inline BBox2f compute_transformed_bbox_fast(Vector2i   const& image_size,
                                              TransformT const& transform_func,
                                              double min_image_size = VW_DEFAULT_MIN_TRANSFORM_IMAGE_SIZE,
                                              double max_image_size = VW_DEFAULT_MAX_TRANSFORM_IMAGE_SIZE) {
    return compute_transformed_bbox_fast(BBox2i(0, 0, image_size[0], image_size[1]), transform_func,
                                                min_image_size, max_image_size);
  }


  // ------------------------
  // class TransformView
  // ------------------------

  /// An image view for transforming an image with an arbitrary mapping functor.
  /// - Basically what this view does is to apply a geometric transform to the
  ///   input i,j coordinates when the accessor operator (i, j) is called.
  template <class ImageT, class TransformT>
  class TransformView : public ImageViewBase<TransformView<ImageT,TransformT> > {

    ImageT     m_image;
    TransformT m_mapper;
    int32      m_width, m_height;

  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef                  pixel_type result_type;
    typedef ProceduralPixelAccessor<TransformView> pixel_accessor;

    // The default constructor creates a tranformed image with the
    // same dimensions as the original.
    TransformView( ImageT const& view, TransformT const& mapper ) :
      m_image(view), m_mapper(mapper), m_width(view.cols()), m_height(view.rows()) {}

    // This constructor allows you to specify the size of the transformed image.
    TransformView( ImageT const& view, TransformT const& mapper, int32 width, int32 height ) :
      m_image(view), m_mapper(mapper), m_width(width), m_height(height) {}

    inline int32 cols  () const { return m_width;          }
    inline int32 rows  () const { return m_height;         }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline result_type operator()( double i, double j, int32 p=0 ) const {
      Vector2 pt_backprojected = m_mapper.reverse(Vector2(i,j));
      return m_image(pt_backprojected[0], pt_backprojected[1], p);
    }

    ImageT     const& child()     const { return m_image;  }
    TransformT const& transform() const { return m_mapper; }

    // \cond INTERNAL
    typedef TransformView<typename ImageT::prerasterize_type, TransformT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i transformed_bbox = m_mapper.reverse_bbox(bbox);
      return prerasterize_type( m_image.prerasterize(transformed_bbox), m_mapper, m_width, m_height );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      if( m_mapper.tolerance() > 0.0 ) {
        ApproximateTransform<TransformT> approx_transform( m_mapper, bbox );
        TransformView<ImageT, ApproximateTransform<TransformT> > approx_view( m_image, approx_transform, m_width, m_height );
        vw::rasterize( approx_view.prerasterize(bbox), dest, bbox );
      }
      else {
        vw::rasterize( prerasterize(bbox), dest, bbox );
      }
    }
    // \endcond
  };

  // \cond INTERNAL
  // Type Traits
  template <class ImplT, class TransformT>
  struct IsFloatingPointIndexable<TransformView<ImplT, TransformT> > : public true_type {};
  // \endcond

  // ------------------------
  // class TransformViewNoData
  // ------------------------

  // An image view for transforming an image with an arbitrary mapping functor and a
  // nodata value.
  template <class ImageT, class TransformT>
  class TransformViewNoData : public ImageViewBase<TransformViewNoData<ImageT,TransformT> > {

    typename ImageT::pixel_type m_nodata_val;
    
    ImageT     m_image;
    TransformT m_mapper;
    int32      m_width, m_height;
    int        m_pixel_buffer;
  public:
    typedef typename ImageT::pixel_type                  pixel_type;
    typedef pixel_type                                   result_type;
    typedef ProceduralPixelAccessor<TransformViewNoData> pixel_accessor;

    // This constructor allows you to specify the size of the
    // transformed image.
    TransformViewNoData( ImageT const& view, TransformT const& mapper,
                         int32 width, int32 height,
                         typename ImageT::pixel_type nodata_val,
                         int pixel_buffer) :
      m_image(view), m_mapper(mapper), m_width(width), m_height(height),
      m_nodata_val(nodata_val), m_pixel_buffer(pixel_buffer) {}

    inline int32 cols  () const { return m_width; }
    inline int32 rows  () const { return m_height; }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline result_type operator()( double i, double j, int32 p=0 ) const {
      Vector2 rv = m_mapper.reverse(Vector2(i,j));
      int b = m_pixel_buffer; // == 2 for bicubic interpolation
      if (rv[0] < b - 1 || rv[0] >= m_image.cols() - b || // out of bounds
          rv[1] < b - 1 || rv[1] >= m_image.rows() - b    // out of bounds
          ) return m_nodata_val;
      return m_image(rv[0], rv[1], p);
    }

    ImageT const& child() const { return m_image; }
    TransformT const& transform() const { return m_mapper; }

    // \cond INTERNAL
    typedef TransformViewNoData<typename ImageT::prerasterize_type,
                                TransformT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i transformed_bbox = m_mapper.reverse_bbox(bbox);
      return prerasterize_type( m_image.prerasterize(transformed_bbox), m_mapper,
                                m_width, m_height, m_nodata_val, m_pixel_buffer);
    }
    template <class DestT> inline void rasterize( DestT const& dest,
                                                  BBox2i const& bbox ) const {
      if( m_mapper.tolerance() > 0.0 ) {
        ApproximateTransform<TransformT> approx_transform( m_mapper, bbox );
        TransformViewNoData<ImageT, ApproximateTransform<TransformT> > approx_view( m_image, approx_transform, m_width, m_height, m_nodata_val, m_pixel_buffer );
        vw::rasterize( approx_view.prerasterize(bbox), dest, bbox );
      }
      else {
        vw::rasterize( prerasterize(bbox), dest, bbox );
      }
    }
    // \endcond
  };

  // \cond INTERNAL
  // Type Traits
  template <class ImplT, class TransformT>
  struct IsFloatingPointIndexable<TransformViewNoData<ImplT, TransformT> > : public true_type {};
  // \endcond

  // -------------------------------------------------------------------------------
  // Functional API
  // -------------------------------------------------------------------------------

  // Apply a Transformation using an arbitary mapping functor.
  //
  // Returns a transformed image view.  The user can specify the type
  // of interpolation and edge extension to be done by supplying the
  // appropriate functors in the last two arguments.  For example:
  //
  // <TT>transform(input, transform_functor, BicubicInterpolation(), ZeroEdgeExtension());</TT>
  //
  // See Interpolation.h and EdgeExtension.h for a list of built-in
  // interpolation and edge extension modes.  It is also possible to
  // write your own.  Again, see the above files for more information.
  template <class ImageT, class TransformT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT> >::type
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func,
                    EdgeT                 const& edge_func,
                    InterpT               const& interp_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
      (interpolate(v, interp_func, edge_func), transform_func);
  }

  // Convenience function: transform with a default interpolation
  // scheme of bilinear interpolation.
  template <class ImageT, class TransformT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func,
                    EdgeT                 const& edge_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
      (interpolate(v, BilinearInterpolation(), edge_func), transform_func);
  }

  // Convenience function: transform with a default scheme of
  // bilinear interpolation and zero edge extension.
  template <class ImageT, class TransformT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func);
  }


  // This variant of transform allows the user to specify the
  // dimensions of the transformed image.  The upper left hand point
  // (0,0) stays fixed.  For a more flexible method of cropping to an
  // arbitrary bounding box, use one of the transform methods defined below.
  template <class ImageT, class TransformT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT const& transform_func,
                    int32 width,
                    int32 height,
                    EdgeT  const& edge_func,
                    InterpT const& interp_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
      (interpolate(v, interp_func, edge_func), transform_func, width, height);
  }

  // Convenience function: transform with a default interpolation
  // scheme of bilinear interpolation. The user can specify the
  // dimensions of the output image.
  template <class ImageT, class TransformT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT const& transform_func,
                    int32 width,
                    int32 height,
                    EdgeT const& edge_func ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
      (interpolate(v, BilinearInterpolation(), edge_func), transform_func, width, height);
  }

  // Convenience function: transform with a default scheme of
  // bilinear interpolation and zero edge extension.  The user can
  // specify the dimensions of the output image.
  template <class ImageT, class TransformT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT const& transform_func,
                    int32 width,
                    int32 height ) {
    return TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
      (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func, width, height);
  }


  // Transform an image and select the bounding box in the
  // transformed space from which to take the new image.  The
  // compute_transformed_bbox() method can be used to compute the
  // bounding box that will fit all of the transformed pixels.
  template <class ImageT, class TransformT, class BBoxRealT, class EdgeT, class InterpT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func,
                    BBox<BBoxRealT,2>     const& bbox,
                    EdgeT                 const& edge_func,
                    InterpT               const& interp_func ) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
                (interpolate(v, interp_func, edge_func), transform_func), bbox);
  }

  // Convenience function: free transform with a default scheme of
  // bilinear interpolation.  The user supplies a bounding box in the
  // transformed space from that determines what pixels will be
  // rasterized.
  template <class ImageT, class TransformT, class BBoxRealT, class EdgeT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func,
                    BBox<BBoxRealT,2>     const& bbox,
                    EdgeT                 const& edge_func ) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TransformT>
                (interpolate(v, BilinearInterpolation(), edge_func), transform_func), bbox);
  }

  // Convenience function: transform with a default scheme of
  // bilinear interpolation and zero edge extension. The user
  // supplies a bounding box in the transformed space from that
  // determines what pixels will be rasterized.
  template <class ImageT, class TransformT, class BBoxRealT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT> >
  inline transform( ImageViewBase<ImageT> const& v,
                    TransformT            const& transform_func,
                    BBox<BBoxRealT,2>     const& bbox ) {
    return crop(TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TransformT>
                (interpolate(v, BilinearInterpolation(), ZeroEdgeExtension()), transform_func), bbox);
  }

  // This variant of transform allows the user to specify the
  // dimensions of the transformed image and a no-data value to be
  // used when the reverse transform returns out-of range pixels.
  template <class ImageT, class TransformT, class EdgeT, class InterpT>
  TransformViewNoData<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
  inline transform_nodata( ImageViewBase<ImageT> const& v,
                           TransformT const& transform_func,
                           int32             width,
                           int32             height,
                           EdgeT      const& edge_func,
                           InterpT    const& interp_func,
                           typename ImageT::pixel_type nodata_val
                           ) {
    return TransformViewNoData<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TransformT>
      (interpolate(v, interp_func, edge_func), transform_func, width, height, nodata_val, interp_func.pixel_buffer);
  }

  // -------------------------------------------------------------------------------
  // Resample
  // -------------------------------------------------------------------------------

  // Resample the image.  The user specifies the scaling factor in x and y.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v,
                   double x_scale_factor,
                   double y_scale_factor,
                   EdgeT   const& edge_func,
                   InterpT const& interp_func ) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor),
                     int(.5+(v.impl().cols()*x_scale_factor)), int(.5+(v.impl().rows()*y_scale_factor)),
                     edge_func, interp_func);
  }

  // Resample the image.  The user specifies the scaling factor in x and y.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, ResampleTransform>
  inline resample( ImageViewBase<ImageT> const& v,
                   double x_scale_factor,
                   double y_scale_factor,
                   EdgeT const& edge_func ) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor),
                     int(.5+(v.impl().cols()*x_scale_factor)), int(.5+(v.impl().rows()*y_scale_factor)),
                     edge_func, vw::BilinearInterpolation());
  }

  // Resample the image.  The user specifies the scaling factor in x and y.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform>
  inline resample( ImageViewBase<ImageT> const& v,
                   double x_scale_factor,
                   double y_scale_factor ) {
    return transform(v, ResampleTransform(x_scale_factor, y_scale_factor),
                     int(.5+(v.impl().cols()*x_scale_factor)), int(.5+(v.impl().rows()*y_scale_factor)),
                     vw::ConstantEdgeExtension(), vw::BilinearInterpolation());
  }

  // Resample the image.  The user specifies the scaling factor in x
  // and y and the dimensions of the output image.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v,
                   double scale_factor,
                   int32 output_width,
                   int32 output_height,
                   EdgeT   const& edge_func,
                   InterpT const& interp_func ) {
    return transform(v, ResampleTransform(scale_factor, scale_factor),
                     output_width, output_height,
                     edge_func, interp_func);
  }

  // Resample the image.  The user specifies the scaling factor in x
  // and y and the dimensions of the output image.
  template <class ImageT, class EdgeT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v,
                   double scale_factor,
                   EdgeT   const& edge_func,
                   InterpT const& interp_func ) {
    return transform(v, ResampleTransform(scale_factor, scale_factor),
                     int(.5+(v.impl().cols()*scale_factor)), int(.5+(v.impl().rows()*scale_factor)),
                     edge_func, interp_func);
  }

  // Resample the image.  The user specifies the scaling factor in x
  // and y and the dimensions of the output image.
  template <class ImageT, class EdgeT>
  typename boost::disable_if<IsScalar<EdgeT>, TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, ResampleTransform> >::type
  inline resample( ImageViewBase<ImageT> const& v,
                   double scale_factor,
                   EdgeT const& edge_func ) {
    return transform(v, ResampleTransform(scale_factor, scale_factor),
                     int(.5+(v.impl().cols()*scale_factor)), int(.5+(v.impl().rows()*scale_factor)),
                     edge_func, vw::BilinearInterpolation());
  }

  // Resample the image.  The user specifies the scaling factor in x
  // and y and the dimensions of the output image.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform>
  inline resample( ImageViewBase<ImageT> const& v,
                   double scale_factor ) {
    return transform(v, ResampleTransform(scale_factor, scale_factor),
                     int(.5+(v.impl().cols()*scale_factor)), int(.5+(v.impl().rows()*scale_factor)),
                     vw::ConstantEdgeExtension(), vw::BilinearInterpolation());
  }


  // -------------------------------------------------------------------------------
  // Resize
  // -------------------------------------------------------------------------------

  // Resize the image.  The user specifies the dimensions of the output image.
  template <class ImageT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, ResampleTransform>
  inline resize( ImageViewBase<ImageT> const& v,
                 int32 output_width,
                 int32 output_height,
                 EdgeT   const& edge_func,
                 InterpT const& interp_func ) {
    return transform(v, ResampleTransform(output_width/(double)v.impl().cols(), output_height/(double)v.impl().rows()),
                     output_width, output_height, edge_func, interp_func );
  }

  // Resize the image.  The user specifies the dimensions of the output image.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, ResampleTransform>
  inline resize( ImageViewBase<ImageT> const& v,
                 int32 output_width,
                 int32 output_height,
                 EdgeT const& edge_func ) {
    return transform( v, ResampleTransform(output_width/(double)v.impl().cols(), output_height/(double)v.impl().rows()),
                      output_width, output_height, edge_func, BilinearInterpolation() );
  }

  // Resize the image.  The user specifies the dimensions of the output image.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform>
  inline resize( ImageViewBase<ImageT> const& v,
                 int32 output_width,
                 int32 output_height ) {
    return transform( v, ResampleTransform(output_width/(double)v.impl().cols(), output_height/(double)v.impl().rows()),
                      output_width, output_height, ConstantEdgeExtension(), BilinearInterpolation() );
  }


  // -------------------------------------------------------------------------------
  // Translate
  // -------------------------------------------------------------------------------

  // Translate the image.  The user specifies the offset in x and y.
  template <class ImageT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, TranslateTransform>
  inline translate( ImageViewBase<ImageT> const& v,
                    double x_offset,
                    double y_offset,
                    EdgeT   const& edge_func,
                    InterpT const& interp_func ) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     edge_func, interp_func);
  }

  // Translate the image.  The user specifies the offset in x and y.
  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, TranslateTransform>
  inline translate( ImageViewBase<ImageT> const& v,
                    double x_offset,
                    double y_offset,
                    EdgeT const& edge_func ) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     edge_func, BilinearInterpolation());
  }

  // Translate the image.  The user specifies the offset in x and y.
  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, TranslateTransform>
  inline translate( ImageViewBase<ImageT> const& v,
                    double x_offset,
                    double y_offset ) {
    return transform(v, TranslateTransform(x_offset, y_offset),
                     v.impl().cols(), v.impl().rows(),
                     ZeroEdgeExtension(), BilinearInterpolation());
  }

  /// Translate the image.  The user specifies the offset in x and y.
  /// This is a special optimized overload for integer offsets.
  template <class ImageT, class EdgeT>
  EdgeExtensionView<ImageT,EdgeT>
  inline translate( ImageViewBase<ImageT> const& im,
                    int32 x_offset,
                    int32 y_offset,
                    EdgeT const& edge_func ) {
    return EdgeExtensionView<ImageT,EdgeT>( im.impl(), -x_offset, -y_offset, im.impl().cols(), im.impl().rows(), edge_func );
  }

  /// Translate the image.  The user specifies the offset in x and y.
  /// This is a special optimized overload for integer offsets.
  template <class ImageT>
  EdgeExtensionView<ImageT,ZeroEdgeExtension>
  inline translate( ImageViewBase<ImageT> const& im,
                    int32 x_offset,
                    int32 y_offset ) {
    return EdgeExtensionView<ImageT,ZeroEdgeExtension>( im.impl(), -x_offset, -y_offset, im.impl().cols(), im.impl().rows(), ZeroEdgeExtension() );
  }


  // -------------------------------------------------------------------------------
  // Rotate
  // -------------------------------------------------------------------------------
  // TODO: Dimensions of destination image.

  // Rotate the image.  The user specifies the angle.
  template <class ImageT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, RotateTransform>
  inline rotate( ImageViewBase<ImageT> const& v,
                 double theta, Vector2 translate,
                 EdgeT   const& edge_func,
                 InterpT const& interp_func ) {
    return transform(v, RotateTransform(theta, translate), /* dims */
                     edge_func, interp_func);
  }

  template <class ImageT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, RotateTransform>
  inline rotate( ImageViewBase<ImageT> const& v,
                 double theta, Vector2 translate,
                 EdgeT const& edge_func ) {
    return transform(v, RotateTransform(theta, translate),
                     edge_func, BilinearInterpolation());
  }

  template <class ImageT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, RotateTransform>
  inline rotate( ImageViewBase<ImageT> const& v,
                 double theta, Vector2 translate=Vector2() ) {
    return transform(v, RotateTransform(theta, translate), /* dims */
                     ZeroEdgeExtension(), BilinearInterpolation());
  }

} // namespace vw

#endif // __VW_IMAGE_TRANSFORM_H__
