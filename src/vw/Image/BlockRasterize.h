// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BlockRasterizeView.h
///
/// Defines a wrapper view that causes an image view to be rasterized
/// one block at a time, rather than all at once, and supports multi-
/// threaded rasterization as well as block-level cacheing.  Even when
/// you don't need those features, in some situations rasterizing one
/// block at a time can dramatically improve performance by reducing
/// memory utilization.
///
#ifndef __VW_IMAGE_BLOCKRASTERIZE_H__
#define __VW_IMAGE_BLOCKRASTERIZE_H__

#include <vw/Core/Cache.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/BlockProcessor.h>

namespace vw {

  /// A wrapper view that rasterizes its child in blocks.
  template <class ImageT>
  class BlockRasterizeView : public ImageViewBase<BlockRasterizeView<ImageT> > {
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::pixel_type result_type;
    typedef ProceduralPixelAccessor<BlockRasterizeView> pixel_accessor;

    BlockRasterizeView( ImageT const& image, Vector2i const& block_size,
                        int num_threads = 0, Cache *cache = NULL )
      : m_child( new ImageT(image) ),
        m_block_size( block_size ),
        m_num_threads( num_threads ),
        m_cache_ptr( cache )
    {
      initialize();
    }

    inline int32 cols() const { return m_child->cols(); }
    inline int32 rows() const { return m_child->rows(); }
    inline int32 planes() const { return m_child->planes(); }
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0, 0 ); }

    inline result_type operator()( int32 x, int32 y, int32 p = 0 ) const {
#if VW_DEBUG_LEVEL > 1
      vw_out(VerboseDebugMessage, "image") << "ImageResourceView rasterizing pixel (" << x << "," << y << "," << p << ")" << std::endl;
#endif
      if ( m_cache_ptr ) {
        // Early-out optimization for single-block resources
        if( m_block_table->size() == 1 ) {
          return ((*m_block_table)[0])->operator()( x, y, p );
        }
        int32 ix = x/m_block_size.x(), iy = y/m_block_size.y();
        return block(ix,iy)->operator()( x-ix*m_block_size.x(), y - iy*m_block_size.y(), p );
      }
      else return (*m_child)(x,y,p);
    }

    ImageT& child() { return *m_child; }
    ImageT const& child() const { return *m_child; }

    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const {
      ImageView<pixel_type> buf( bbox.width(), bbox.height(), planes() );
      rasterize( buf, bbox );
      return CropView<ImageView<pixel_type> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),cols(),rows()) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      RasterizeFunctor<DestT> rasterizer( *this, dest, bbox.min() );
      BlockProcessor<RasterizeFunctor<DestT> > process( rasterizer, m_block_size, m_num_threads );
      process(bbox);
    }

  private:
    // These function objects are spawned to rasterize the child image.
    // One functor is created per child thread, and they are called
    // in succession with bounding boxes that are each contained
    // within one block.
    template <class DestT>
    class RasterizeFunctor {
      BlockRasterizeView const& m_view;
      DestT const& m_dest;
      Vector2i m_offset;
    public:
      RasterizeFunctor( BlockRasterizeView const& view, DestT const& dest, Vector2i const& offset )
        : m_view(view), m_dest(dest), m_offset(offset) {}
      void operator()( BBox2i const& bbox ) const {
#if VW_DEBUG_LEVEL > 1
        vw_out(VerboseDebugMessage, "image") << "BlockRasterizeView::RasterizeFunctor( " << bbox << " )" << std::endl;
#endif
        if( m_view.m_cache_ptr ) {
          int32 ix=bbox.min().x()/m_view.m_block_size.x(), iy=bbox.min().y()/m_view.m_block_size.y();
#if VW_DEBUG_LEVEL > 1
          int32 maxix=(bbox.max().x()-1)/m_view.m_block_size.x(), maxiy=(bbox.max().y()-1)/m_view.m_block_size.y();
          if( maxix != ix || maxiy != iy ) {
            vw_throw(LogicErr() << "BlockRasterizeView::RasterizeFunctor: bbox spans more than one cache block!");
          }
#endif
          m_view.block(ix,iy)->rasterize( crop( m_dest, bbox-m_offset ), bbox-Vector2i(ix*m_view.m_block_size.x(),iy*m_view.m_block_size.y()) );
        }
        else m_view.child().rasterize( crop( m_dest, bbox-m_offset ), bbox );
      }
    };

    // Allows RasterizeFunctor to access cache-related members.
    template <class DestT> friend class RasterizeFunctor;

    // These objects rasterize a full block of image data to be
    // stored in the cache.
    class BlockGenerator {
      boost::shared_ptr<ImageT> m_child;
      BBox2i m_bbox;
    public:
      typedef ImageView<pixel_type> value_type;

      BlockGenerator( boost::shared_ptr<ImageT> const& child, BBox2i const& bbox )
        : m_child( child ), m_bbox( bbox ) {}

      size_t size() const {
        return m_bbox.width() * m_bbox.height() * m_child->planes() * sizeof(pixel_type);
      }

      boost::shared_ptr<ImageView<pixel_type> > generate() const {
        boost::shared_ptr<ImageView<pixel_type> > ptr( new ImageView<pixel_type>( m_bbox.width(), m_bbox.height(), m_child->planes() ) );
        m_child->rasterize( *ptr, m_bbox );
        return ptr;
      }
    };

    void initialize() {
      if( m_block_size.x() <= 0 || m_block_size.y() <= 0 ) {
        const int32 default_blocksize = 2*1024*1024; // 2 megabytes
        // XXX Should the default block configuration be different for
        // very wide images?  Either way we will guess wrong some of
        // the time, so advanced users will have to know what they're
        // doing in any case.
        int32 block_rows = default_blocksize / (planes()*cols()*int32(sizeof(pixel_type)));
        if( block_rows < 1 ) block_rows = 1;
        else if( block_rows > rows() ) block_rows = rows();
        m_block_size = Vector2i( cols(), block_rows );
      }
      if( m_cache_ptr ) {
        m_table_width = (cols()-1) / m_block_size.x() + 1;
        m_table_height = (rows()-1) / m_block_size.y() + 1;
        m_block_table.reset( new std::vector<Cache::Handle<BlockGenerator> >( m_table_width * m_table_height ) );
        BBox2i view_bbox(0,0,cols(),rows());
        for( int32 iy=0; iy<m_table_height; ++iy ) {
          for( int32 ix=0; ix<m_table_width; ++ix ) {
            BBox2i bbox( ix*m_block_size.x(), iy*m_block_size.y(), m_block_size.x(), m_block_size.y() );
            bbox.crop( view_bbox );
            block(ix,iy) = m_cache_ptr->insert( BlockGenerator( m_child, bbox ) );
          }
        }
      }
    }

    Cache::Handle<BlockGenerator>& block( int ix, int iy ) const {
      if( ix<0 || ix>=m_table_width || iy<0 || iy>=m_table_height )
        vw_throw( ArgumentErr() << "BlockRasterizeView: Block indices out of bounds, (" << ix
                  << "," << iy << ") of (" << m_table_width << "," << m_table_height << ")" );
      return (*m_block_table)[ix+iy*m_table_width];
    }

    // We store this by shared pointer so it doesn't move when we copy
    // the BlockRasterizeView, since the BlockGenerators point to it.
    boost::shared_ptr<ImageT> m_child;
    Vector2i m_block_size;
    int32 m_num_threads;
    Cache *m_cache_ptr;
    int m_table_width, m_table_height;
    // We store this by shared pointer so copying a BlockRasterizeView
    // (i.e. to promote its scope) is not as expensive an operation.
    boost::shared_ptr<std::vector<Cache::Handle<BlockGenerator> > > m_block_table;
  };

  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_rasterize( ImageViewBase<ImageT> const& image, Vector2i const& block_size, int num_threads = 0 ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads );
  }

  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_cache( ImageViewBase<ImageT> const& image, Vector2i const& block_size, int num_threads = 0 ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads, &vw_system_cache() );
  }

  template <class ImageT>
  inline BlockRasterizeView<ImageT> block_cache( ImageViewBase<ImageT> const& image, Vector2i const& block_size, int num_threads, Cache& cache ) {
    return BlockRasterizeView<ImageT>( image.impl(), block_size, num_threads, &cache );
  }

} // namespace vw

#endif // __VW_IMAGE_BLOCKRASTERIZE_H__
