// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ViewImageResource.h
///
/// This ImageResource can wrap any vision workbench image view so
/// that it can be presented as an ImageResource to other
/// subsystems.
///
#ifndef __VW_IMAGE_VIEW_IMAGERESOURCE_H__
#define __VW_IMAGE_VIEW_IMAGERESOURCE_H__

#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

namespace vw {

  /// \cond INTERNAL
  // Base class definition
  class ViewImageResourceBase {
  public:
    virtual ~ViewImageResourceBase() {}

    /// Returns the number of columns in an image resource.
    virtual int32 cols() const = 0;

    /// Returns the number of rows in an image resource.
    virtual int32 rows() const = 0;

    /// Returns the number of planes in an image resource.
    virtual int32 planes() const = 0;

    /// Returns the number of channels in a image resource.
    virtual int32 channels() const = 0;

    /// Returns the native pixel format of the resource.
    virtual PixelFormatEnum pixel_format() const = 0;

    /// Returns the native channel type of the resource.
    virtual ChannelTypeEnum channel_type() const = 0;

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const = 0;

    /// Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox ) = 0;

    virtual char* data() const = 0;
  };



  // This set of classes helps us to extract the direct data access
  // member (if available) from the view.
  template <class ViewT> struct ViewDataAccessor {
    static char* data(ViewT const& /*view*/) {
      vw_throw(NoImplErr() << "ViewDataAccessor data() failed. This view does not support direct data access.");
      return NULL; // never reached
    }
  };

  // Currently, the ImageView<> class is the only one that supports
  // direct access.
  template<class PixelT> struct ViewDataAccessor<ImageView<PixelT> > {
    static char* data(ImageView<PixelT> const& view) { return (char*)(view.data()); }
  };


  // ViewImageResourceImpl class implementation
  template <class ViewT>
  class ViewImageResourceImpl : public ViewImageResourceBase {
  private:
    ViewT m_view;

  public:
    ViewImageResourceImpl( ImageViewBase<ViewT> const& view ) :
      m_view(view.impl()) {}

    virtual ~ViewImageResourceImpl() {}

    /// Returns the number of columns in an image resource.
    virtual int32 cols() const { return m_view.cols(); }

    /// Returns the number of rows in an image resource.
    virtual int32 rows() const { return m_view.rows(); }

    /// Returns the number of planes in an image resource.
    virtual int32 planes() const { return m_view.planes(); }

    /// Returns the number of planes in an image resource.
    virtual int32 channels() const { return m_view.channels(); }

    /// Returns the native pixel format of the resource.
    virtual PixelFormatEnum pixel_format() const { return m_view.pixel_format(); }

    /// Returns the native channel type of the resource.
    virtual ChannelTypeEnum channel_type() const { return m_view.channel_type(); }

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const {
      VW_ASSERT( bbox.min().x() >= 0 && bbox.min().y() >= 0 &&
                 bbox.max().x() <= cols() && bbox.max().y() <= rows(),
                 ArgumentErr() << "ViewImageResource::read(): bbox exceeds view dimensions." );
      ImageView<typename ViewT::pixel_type> region = crop(m_view, bbox);
      ImageBuffer src = region.buffer();
      convert(buf, src);
    }

    /// Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& /*buf*/, BBox2i const& /*bbox*/ ) {
      vw_throw(NoImplErr() << "ViewImageResource::write() failed.  ViewImageResource is read-only.");
    }

    virtual char* data() const {
      return ViewDataAccessor<ViewT>::data(m_view);
    }
  };

  /// Base class from which specific image resources derive.
  class ViewImageResource : public ImageResource {
    boost::shared_ptr<ViewImageResourceBase> m_rsrc;
    Vector2i m_block_size;

  public:

    template <class ViewT>
    explicit ViewImageResource(ImageViewBase<ViewT> const& view, Vector2i block_size) :
      m_rsrc( new ViewImageResourceImpl<ViewT>(view) ), m_block_size(block_size) {}

    template <class ViewT>
    explicit ViewImageResource(ImageViewBase<ViewT> const& view) :
      m_rsrc( new ViewImageResourceImpl<ViewT>(view) ),
      m_block_size(Vector2i(view.impl().cols(), view.impl().rows())) {}

    virtual ~ViewImageResource() {};

    /// Returns the number of columns in an image resource.
    virtual int32 cols() const { return m_rsrc->cols(); }

    /// Returns the number of rows in an image resource.
    virtual int32 rows() const { return m_rsrc->rows(); }

    /// Returns the number of planes in an image resource.
    virtual int32 planes() const { return m_rsrc->planes(); }

    /// Returns the number of channels in a image resource.
    int32 channels() const { return num_channels( pixel_format() ); }

    /// Returns the native pixel format of the resource.
    virtual PixelFormatEnum pixel_format() const { return m_rsrc->pixel_format(); }

    /// Returns the native channel type of the resource.
    virtual ChannelTypeEnum channel_type() const { return m_rsrc->channel_type(); }

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const {
      return m_rsrc->read(buf, bbox);
    }

    /// Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox ) {
      return m_rsrc->write(buf, bbox);
    }

    /// Returns the optimal block size/alignment for partial reads or writes.
    virtual Vector2i block_size() const { return m_block_size; }

    /// Set the preferred block size/alignment for partial reads or writes.
    virtual void set_block_size( Vector2i const& size ) { m_block_size = size; }

    /// Force any changes to be written to the resource.
    virtual void flush() {}

    char* data() const { return m_rsrc->data(); }
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGERESOURCE_H__
