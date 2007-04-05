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

/// \file ImageResourceView.h
///
/// A block-cached read-only view of a general image resource.
///
#ifndef __VW_IMAGE_IMAGERESOURCEVIEW_H__
#define __VW_IMAGE_IMAGERESOURCEVIEW_H__

#include <string>
#include <map>

#include <vw/Core/Cache.h>
#include <vw/Core/Debugging.h>
#include <vw/Image/ImageView.h>

namespace vw {

  /// A view of an image resourcs.
  template <class PixelT>
  class ImageResourceView : public ImageViewBase<ImageResourceView<PixelT> >
  {

    // BlockGenerator
    class BlockGenerator {
    public:
      boost::shared_ptr<ImageResource> m_res_ptr;
      BBox2i m_bbox;
    public:
      typedef ImageView<PixelT> value_type;

      BlockGenerator( boost::shared_ptr<ImageResource> res_ptr, BBox2i bbox )
        : m_res_ptr( res_ptr ), m_bbox( bbox ) {}
      
      size_t size() const {
        return m_bbox.width() * m_bbox.height() * m_res_ptr->planes() * sizeof(PixelT);
      }
      
      boost::shared_ptr<ImageView<PixelT> > generate() const {

        // If the user has requested a single-plane, compound pixel
        // type, but the file is a multi-plane, scalr pixel file
        // without any pixel semantics, we force the resource to a
        // single plane so that convert() can convert from a
        // multiplane file to compound pixel type image.
        int planes = m_res_ptr->planes();
        if (PixelNumChannels<PixelT>::value != 1)
          planes = 1;
        
        boost::shared_ptr<ImageView<PixelT> > ptr( new ImageView<PixelT>( m_bbox.width(), m_bbox.height(), planes ) );
        m_res_ptr->read( ptr->buffer(), m_bbox );
        return ptr;
      }
    };

    boost::shared_ptr<ImageResource> r;

    Cache& m_cache;
    bool m_enable_cache;
    Vector2i m_block_size;
    typedef std::map<std::pair<int,int>, Cache::Handle<BlockGenerator> > block_table_type;
    mutable block_table_type m_block_table;

    void initialize() {
      if( unsigned(m_block_size.x())==cols() && m_block_size.y()==1 ) {
        // Group scanlines into 16K chunks for efficiency
        const size_t blocksize = 16384;
        if( cols()*rows()*planes()*sizeof(PixelT) < blocksize ) {
          m_block_size.y()=rows();
        }
        else {
          size_t line_size = cols()*planes()*sizeof(PixelT);
          m_block_size.y() = blocksize/line_size;
          if( m_block_size.y() == 0 ) m_block_size.y() = 1;
        }
      }
      int maxx=(cols()-1)/m_block_size.x();
      int maxy=(rows()-1)/m_block_size.y();
      BBox2i view_bbox(0,0,cols(),rows());
      for( int y=0; y<=maxy; ++y ) {
        for( int x=0; x<=maxx; ++x ) {
          BBox2i bbox( x*m_block_size.x(), y*m_block_size.y(), m_block_size.x(), m_block_size.y() );
          bbox.crop( view_bbox );
          Cache::Handle<BlockGenerator> handle = m_cache.insert( BlockGenerator( r, bbox ) );
          m_block_table[std::make_pair(x,y)] = handle;
        }
      }
    }

  public:
    /// The pixel type of the view.
    typedef PixelT pixel_type;
    typedef PixelT const& result_type;

    /// The view's pixel accessor type.
    typedef ProceduralPixelAccessor<ImageResourceView<PixelT> > pixel_accessor;
    
    /// Constructs an ImageResourceView of the given resource.  Takes
    /// ownership of the resource object (i.e. deletes it when it's
    /// done using it).
    ImageResourceView( ImageResource *resource, bool cache=true )
      : r( resource ),
        m_cache(Cache::system_cache()),
        m_enable_cache(cache),
        m_block_size( r->native_block_size() )
    { if(cache) initialize(); }

    /// Constructs an ImageResourceView of the given resource using the
    /// specified cache area.  Takes ownership of the resource object
    /// (i.e. deletes it when it's done using it).
    ImageResourceView( ImageResource *resource, Cache& cache )
      : r( resource ),
        m_cache(cache),
        m_enable_cache(true),
        m_block_size( r->native_block_size() )
    { initialize(); }

    ~ImageResourceView() {}
    
    /// Returns the number of columns in the image.
    inline unsigned cols() const { return r->cols(); }

    /// Returns the number of rows in the image.
    inline unsigned rows() const { return r->rows(); }

    /// Returns the number of planes in the image.
    inline unsigned planes() const { return 1; }
    
    /// Returns the pixel at the given position in the given plane.
    result_type operator()( unsigned x, unsigned y, unsigned /*plane*/=1 ) const {
      if( ! m_enable_cache )
        vw_throw( LogicErr() << "Non-cacheing ImageResourceViews do not support per-pixel access" );
      int ix = x/m_block_size.x(), iy = y/m_block_size.y();
      return m_block_table[std::make_pair(ix,iy)]->operator()( x-ix*m_block_size.x(), y - iy*m_block_size.y() );
    }
    
    /// Returns a pixel_accessor pointing to the origin.
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    
    ImageResource const* resource() const { return r.get(); }

    /// \cond INTERNAL
    typedef ImageResourceView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      if( !m_enable_cache ) return read_image( dest, *r, bbox );
      int ix0=bbox.min().x()/m_block_size.x(), iy=bbox.min().y()/m_block_size.y();
      int maxix=(bbox.max().x()-1)/m_block_size.x(), maxiy=(bbox.max().y()-1)/m_block_size.y();
      for( ; iy <= maxiy; ++iy ) {
        for( int ix = ix0; ix <= maxix; ++ix ) {
          BBox2i block_bbox( ix*m_block_size.x(), iy*m_block_size.y(), m_block_size.x(), m_block_size.y() );
          block_bbox.crop( bbox );
          m_block_table[std::make_pair(ix,iy)]->rasterize( crop( dest, block_bbox-bbox.min() ),
                                                           block_bbox-Vector2i(ix*m_block_size.x(),iy*m_block_size.y()) );
        }
      }
    }
    /// \endcond
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGERESOURCEVIEW_H__
