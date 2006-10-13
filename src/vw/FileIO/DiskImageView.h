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

/// \file DiskImageView.h
///
/// This is a barely-functional first cut at on-disk image 
/// support, currently used by the panorama blending code.
/// It is read-only, supports only VIL-based file types, 
/// doesn't support multi-plane images, doesn't support 
/// random pixel access, etc.  Buyer beware.
///
#ifndef __VW_FILEIO_DISK_IMAGE_VIEW_H__
#define __VW_FILEIO_DISK_IMAGE_VIEW_H__

#include <string>
#include <map>

#include <vw/Core/Cache.h>
#include <vw/Core/Debugging.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  /// A view of an image on disk.
  template <class PixelT>
  class DiskImageView : public ImageViewBase<DiskImageView<PixelT> >
  {

    // BlockGenerator
    class BlockGenerator {
    public:
      boost::shared_ptr<DiskImageResource> m_res_ptr;
      BBox2i m_bbox;
    public:
      typedef ImageView<PixelT> value_type;

      BlockGenerator( boost::shared_ptr<DiskImageResource> res_ptr, BBox2i bbox )
        : m_res_ptr( res_ptr ), m_bbox( bbox ) {}
      
      size_t size() const {
        return m_bbox.width() * m_bbox.height() * m_res_ptr->planes() * sizeof(PixelT);
      }
      
      boost::shared_ptr<ImageView<PixelT> > generate() const {
        vw_out(DebugMessage) << "DiskImageView reading block " << m_bbox << " from " << m_res_ptr->filename() << std::endl;
        boost::shared_ptr<ImageView<PixelT> > ptr( new ImageView<PixelT>( m_bbox.width(), m_bbox.height(), m_res_ptr->planes() ) );
        m_res_ptr->read_generic( GenericImageBuffer(*ptr), m_bbox );
        return ptr;
      }
    };

    boost::shared_ptr<DiskImageResource> r;

    Cache& m_cache;
    bool m_enable_cache;
    Vector2i m_block_size;
    typedef std::map<std::pair<int,int>, Cache::Handle<BlockGenerator> > block_table_type;
    mutable block_table_type m_block_table;

    void initialize() {
      if( m_block_size.x()==cols() && m_block_size.y()==1 ) {
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
    typedef ProceduralPixelAccessor<DiskImageView<PixelT> > pixel_accessor;
    
    /// Constructs a DiskImageView of the given file on disk.
    DiskImageView( std::string const& filename, bool cache=true )
      : r( DiskImageResource::open( filename ) ),
        m_cache(Cache::system_cache()),
        m_enable_cache(cache),
        m_block_size( r->native_read_block_size() )
    { if(cache) initialize(); }

    /// Constructs a DiskImageView of the given file on disk 
    /// using the specified cache area.
    DiskImageView( std::string const& filename, Cache& cache )
      : r( DiskImageResource::open( filename ) ),
        m_cache(cache),
        m_enable_cache(true),
        m_block_size( r->native_read_block_size() )
    { initialize(); }

    ~DiskImageView() {}
    
    /// Returns the number of columns in the image.
    inline unsigned cols() const { return r->cols(); }

    /// Returns the number of rows in the image.
    inline unsigned rows() const { return r->rows(); }

    /// Returns the number of planes in the image.
    inline unsigned planes() const { return 1; }
    
    /// Returns the pixel at the given position in the given plane.
    result_type operator()( unsigned i, unsigned j, unsigned plane=1 ) const {
      throw NoImplErr() << "DiskImageView does not currently support random pixel access";
    }
    
    /// Returns a pixel_accessor pointing to the origin.
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    
    std::string filename() const { return r->filename(); }

    /// \cond INTERNAL
    typedef DiskImageView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      if( !m_enable_cache ) return r->read( dest, bbox );
      int ix0=bbox.min().x()/m_block_size.x(), iy=bbox.min().y()/m_block_size.y();
      int maxix=(bbox.max().x()-1)/m_block_size.x(), maxiy=(bbox.max().y()-1)/m_block_size.y();
      vw_out(VerboseDebugMessage) << "DiskImageView rasterizing region " << bbox << " (blocks (" << 
        ix0 << "," << iy << ")-(" << maxix << "," << maxiy << ")) from " << filename() << std::endl;
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

#endif // __VW_FILEIO_DISK_IMAGE_VIEW_H__
