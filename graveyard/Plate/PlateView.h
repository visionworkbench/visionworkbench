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
    boost::shared_ptr<ReadOnlyPlateFile> m_platefile;
    int m_current_level;
    TransactionOrNeg m_transaction_id;

  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ProceduralPixelAccessor<PlateView> pixel_accessor;

    PlateView(const Url& url)
      : m_platefile( new PlateFile(url) ),
        m_current_level(m_platefile->num_levels()-1),
        m_transaction_id(-1)
    { }

    PlateView(boost::shared_ptr<PlateFile> plate)
      : m_platefile( plate ),
        m_current_level(m_platefile->num_levels()-1),
        m_transaction_id(-1)
    { }

    // Standard ImageView interface methods
    int32 cols() const {
      return (1 << m_current_level) * m_platefile->default_tile_size();
    }

    int32 rows() const {
      return (1 << m_current_level) * m_platefile->default_tile_size();
    }

    int32 planes() const { return 1; }

    void set_level(int level) {
      if (level < 0 || level >= num_levels())
        vw_throw(ArgumentErr() << "PlateView::set_level() -- invalid level");
      m_current_level = level;
    }

    void set_transaction(TransactionOrNeg t) {
      m_transaction_id = t;
    }

    int num_levels() const { return m_platefile->num_levels(); }

    std::list<TileHeader>
    search_for_tiles( BBox2i image_bbox ) const {
      const float tile_size = m_platefile->default_tile_size();
      BBox2i query_region;
      query_region.min() = floor(Vector2f(image_bbox.min())/tile_size);
      query_region.max() = ceil(Vector2f(image_bbox.max())/tile_size);
      return m_platefile->search_by_region(m_current_level, query_region, TransactionRange(m_transaction_id));
    }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(int32 /*x*/, int32 /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "PlateView::operator() -- not yet implemented.");
    }

    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    //typedef CropView<TransformView<InterpolationView<EdgeExtensionView<CropView<ImageView<pixel_type> >, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {
      const int32 tile_size = m_platefile->default_tile_size();

      std::list<TileHeader> tileheaders =
        search_for_tiles( bbox );

      if ( tileheaders.empty() )
        return crop( ImageView<pixel_type>(constant_view(pixel_type(),
                                                         bbox.width(), bbox.height())),
                     BBox2i(-bbox.min().x(), -bbox.min().y(), this->cols(), this->rows()) );

      // Create an image of the appropriate size to rasterize tiles into.
      ImageView<pixel_type> level_image(bbox.width(),bbox.height());

      // Access the tiles needed for this level and copy them into place
      BOOST_FOREACH( TileHeader const& theader, tileheaders ) {
        ImageView<PixelT> tile;
        m_platefile->read( tile, theader.col(), theader.row(), theader.level(),
                           theader.transaction_id(), true );

        BBox2i src_bbox_cropped( tile_size*theader.col(), tile_size*theader.row(),
                                    tile_size, tile_size );
        src_bbox_cropped.crop( bbox );

        BBox2i dst_bbox_cropped = src_bbox_cropped;

        src_bbox_cropped.min() -= Vector2i(tile_size*theader.col(),
                                           tile_size*theader.row());
        src_bbox_cropped.max() -= Vector2i(tile_size*theader.col(),
                                           tile_size*theader.row());
        dst_bbox_cropped.min() -= bbox.min();
        dst_bbox_cropped.max() -= bbox.min();

        crop( level_image, dst_bbox_cropped ) =
          crop( tile, src_bbox_cropped );
      }

      return crop( level_image,
                   BBox2i(-bbox.min().x(), -bbox.min().y(), this->cols(), this->rows()) );

    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_PLATEVIEW__
