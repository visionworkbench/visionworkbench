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

  // This set of classes helps us to extract the direct data access
  // member (if available) from the view.
  template <class ViewT> struct ViewDataAccessor {
    static boost::shared_array<const uint8> data(ViewT const& /*view*/) {
      vw_throw(NoImplErr() << "ViewDataAccessor native_ptr() failed. This view does not support direct data access.");
      return NULL; // never reached
    }
  };

  // Currently, the ImageView<> class is the only one that supports
  // direct access.
  template<class PixelT> struct ViewDataAccessor<ImageView<PixelT> > {
    static boost::shared_array<const uint8> data(ImageView<PixelT> const& view) {
      return boost::shared_array<const uint8>(reinterpret_cast<const uint8*>(view.data()), NOP());
    }
  };


  template <class ViewT>
  class ViewImageResourceImpl : public SrcImageResource {
  private:
    ViewT m_view;

  public:
    ViewImageResourceImpl( ImageViewBase<ViewT> const& view ) :
      m_view(view.impl()) {}

    virtual ~ViewImageResourceImpl() {}

    virtual ImageFormat format() const { return m_view.format(); }

    virtual bool has_block_write() const  {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read() const   {return false;}
    virtual bool has_nodata_read() const  {return false;}

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const {
      VW_ASSERT( bbox.min().x() >= 0 && bbox.min().y() >= 0 &&
                 bbox.max().x() <= cols() && bbox.max().y() <= rows(),
                 ArgumentErr() << "ViewImageResource::read(): bbox exceeds view dimensions." );
      ImageView<typename ViewT::pixel_type> region = crop(m_view, bbox);
      ImageBuffer src = region.buffer();
      convert(buf, src);
    }

    virtual boost::shared_array<const uint8> native_ptr() const {
      return ViewDataAccessor<ViewT>::data(m_view);
    }
  };

  /// Base class from which specific image resources derive.
  class ViewImageResource : public SrcImageResource {
    boost::shared_ptr<SrcImageResource> m_rsrc;
    Vector2i m_block_size;

  public:

    template <class ViewT>
    explicit ViewImageResource(ImageViewBase<ViewT> const& view, Vector2i block_size) :
      m_rsrc( new ViewImageResourceImpl<ViewT>(view) ), m_block_size(block_size) {}

    template <class ViewT>
    explicit ViewImageResource(ImageViewBase<ViewT> const& view) :
      m_rsrc( new ViewImageResourceImpl<ViewT>(view) ),
      m_block_size(Vector2i(view.impl().cols(), view.impl().rows())) {}

    virtual ImageFormat format() const { return m_rsrc->format(); }

    /// Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const {
      return m_rsrc->read(buf, bbox);
    }

    virtual bool has_block_write() const  {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read() const   {return true;}
    virtual bool has_nodata_read() const  {return false;}

    /// Returns the optimal block size/alignment for partial reads
    virtual Vector2i block_read_size() const { return m_block_size; }

    /// Set the preferred block size/alignment for partial reads
    virtual void set_block_read_size( Vector2i const& size ) { m_block_size = size; }

    /// Force any changes to be written to the resource.
    virtual void flush() {}

    boost::shared_array<const uint8> native_ptr() const { return m_rsrc->native_ptr(); }
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGERESOURCE_H__
