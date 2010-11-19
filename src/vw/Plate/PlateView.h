// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATEVIEW__
#define __VW_PLATE_PLATEVIEW__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Image/Transform.h>
#include <boost/foreach.hpp>

namespace vw {
namespace platefile {

  /// An image view for accessing tiles from a plate file.  Tiles are
  /// cached by this view to increase read speeds.
  template <class PixelT>
  class PlateView : public ImageViewBase<PlateView<PixelT> > {
    boost::shared_ptr<PlateFile> m_platefile;
    int m_current_level;

  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ProceduralPixelAccessor<PlateView> pixel_accessor;

    PlateView(const Url& url)
      : m_platefile( new PlateFile(url) ),
        m_current_level(m_platefile->num_levels()-1)
    { }

    PlateView(boost::shared_ptr<PlateFile> plate)
      : m_platefile( plate ),
        m_current_level(m_platefile->num_levels()-1)
    { }

    // Standard ImageView interface methods
    int32 cols() const {
      return (1 << (m_platefile->num_levels()-1)) * m_platefile->default_tile_size();
    }

    int32 rows() const {
      return (1 << (m_platefile->num_levels()-1)) * m_platefile->default_tile_size();
    }

    int32 planes() const { return 1; }

    void set_level(int level) {
      if (level < 0 || level >= num_levels())
        vw_throw(ArgumentErr() << "PlateView::set_level() -- invalid level");
      m_current_level = level;
    }

    int num_levels() const {
      return m_platefile->num_levels();
    }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(int32 /*x*/, int32 /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "PlateView::operator() -- not yet implemented.");
    }

    // Returns the number of tiles in a certain image square
    std::list<TileHeader>
    search_for_tiles( BBox2i const& image_bbox, BBox2i& level_bbox,
                      BBox2i& aligned_level_bbox ) const {
      // Compute the bounding box at the current level.
      int level_difference = (m_platefile->num_levels()-1) - m_current_level;
      level_bbox = image_bbox / pow(2,level_difference);

      // Grow that bounding box to align with tile boundaries
      aligned_level_bbox = level_bbox;
      aligned_level_bbox.min() = ( (level_bbox.min() / m_platefile->default_tile_size())
                                  * m_platefile->default_tile_size() );
      aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) /
                                                  float(m_platefile->default_tile_size()) ))
                                       * m_platefile->default_tile_size() );
      aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) /
                                                  float(m_platefile->default_tile_size()) ))
                                       * m_platefile->default_tile_size() );

      // Searching for tile headers
      return m_platefile->search_by_region( m_current_level,
                                            aligned_level_bbox/m_platefile->default_tile_size(),
                                            -1, -1, 0 );
    }

    // Simplified access
    std::list<TileHeader>
    search_for_tiles( BBox2i const& image_bbox ) const {
      BBox2i level_bbox, aligned_level_bbox;
      return search_for_tiles( image_bbox, level_bbox,
                               aligned_level_bbox );
    }

    /// \cond INTERNAL
    typedef CropView<TransformView<InterpolationView<EdgeExtensionView<CropView<ImageView<pixel_type> >, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {

      // Compute the bounding box at the current level.
      int level_difference = (m_platefile->num_levels()-1) - m_current_level;
      BBox2i level_bbox, aligned_level_bbox;

      std::list<TileHeader> tileheaders =
        search_for_tiles( bbox, level_bbox, aligned_level_bbox );

      // Create an image of the appropriate size to rasterize tiles into.
      ImageView<pixel_type> level_image(aligned_level_bbox.width(), aligned_level_bbox.height());

      // Access the tiles needed for this level and copy them into place
      BOOST_FOREACH( TileHeader const& theader, tileheaders ) {
        ImageView<PixelT> tile;
        m_platefile->read( tile, theader.col(), theader.row(), theader.level(),
                           theader.transaction_id(), true );
        BBox2i tile_bbox( m_platefile->default_tile_size()*theader.col()-aligned_level_bbox.min().x(),
                          m_platefile->default_tile_size()*theader.row()-aligned_level_bbox.min().y(),
                          m_platefile->default_tile_size(),
                          m_platefile->default_tile_size() );
        crop( level_image, tile_bbox ) = tile;
      }

      // Crop the output to the original requested bbox, resample it, and return.
      BBox2i output_bbox( level_bbox.min().x() - aligned_level_bbox.min().x(),
                          level_bbox.min().y() - aligned_level_bbox.min().y(),
                          level_bbox.width(), level_bbox.height() );

      return crop(transform( crop(level_image, output_bbox),
                             ResampleTransform( pow(2, level_difference),
                                                pow(2, level_difference) ),
                             ConstantEdgeExtension(), BilinearInterpolation() ),
                  BBox2i(-bbox.min().x(), -bbox.min().y(), this->cols(), this->rows()));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_PLATEVIEW__
