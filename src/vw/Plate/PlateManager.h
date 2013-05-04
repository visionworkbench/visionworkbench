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


#ifndef __VW_PLATE_PLATEMANAGER_H__
#define __VW_PLATE_PLATEMANAGER_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Filter.h>
#include <boost/foreach.hpp>

namespace vw {

  template <typename PixelT> class ImageViewRef;

namespace platefile {

  namespace detail {
    class RememberCallback;
  }

  template <class ViewT>
  class WritePlateFileTask;

  // The Tile Entry is used to keep track of the bounding box of
  // tiles and their location in the grid.
  struct TileInfo {
    int i, j;
    BBox2i bbox;
    TileInfo(int i, int j, BBox2i const& bbox) : i(i), j(j), bbox(bbox) {}
  };

  // Functor for clamping alpha to channel range
  struct ClampAlpha : UnaryReturnSameType {
    template <class T>
    T operator()( T const& px ) const {
      typedef typename CompoundChannelType<T>::type ChT;
      typedef ChannelRange<ChT> ChR;
      T copy = px;
      if ( alpha_channel( copy ) > ChR::max() )
        alpha_channel( copy ) = ChR::max();
      else if ( alpha_channel( copy ) < ChR::min() )
        alpha_channel( copy ) = ChR::min();
      return copy;
    }
  };

  template <class PixelT>
  class PlateManager {
  protected:
    boost::shared_ptr<PlateFile> m_platefile;

    virtual void affected_tiles( BBox2i const& image_size,
                                 TransformRef const& tx, int tile_size,
                                 int level, std::list<TileInfo>& tiles ) const;

    virtual void transform_image( cartography::GeoReference const& georef,
                                  ImageViewRef<PixelT>& image,
                                  TransformRef& txref, int& level ) const = 0;

    void slow_mipmap( uint32 level, std::list<TileHeader>& src_hdrs, bool preblur, const detail::RememberCallback& pc) const;
    void fast_mipmap( uint32 starting_level, int32 stopping_level, std::list<TileHeader>& src_hdrs, bool preblur, const detail::RememberCallback& pc) const;
    uint64 calc_cache_tile_count() const;

  public:
    PlateManager(boost::shared_ptr<PlateFile> platefile) : m_platefile(platefile) {}

    virtual ~PlateManager() {}

    // mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
    //
    // starting_level -- select the pyramid level from which to carry
    //                   out mipmapping
    // bbox -- bounding box (in terms of tiles) containing the tiles that need
    //         to be mipmapped at starting_level. Use to specify affected
    //         tiles.
    // input_transaction_id -- to use when reading
    // output_transaction_id -- to use when writing
    virtual void mipmap(uint32 starting_level, BBox2i const& bbox, TransactionOrNeg input_transaction_id, bool preblur, const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(), uint32 stopping_level=0) const;

    // Provides user a georeference for a particular level of the pyramid
    virtual cartography::GeoReference georeference( int level ) const = 0;

    // Provide the user a platemanager of type mode
    static PlateManager<PixelT>* make( std::string const& mode,
                                       boost::shared_ptr<PlateFile> platefile );

    // Adds an image to the plate file.
    template <class ViewT>
    void insert( ImageViewBase<ViewT> const& imagebase,
                 std::string const& description,
                 TransactionOrNeg transaction_id_override,
                 cartography::GeoReference const& input_georef,
                 bool tweak_settings_for_terrain, bool /*verbose*/ = false,
                 const ProgressCallback &progress = ProgressCallback::dummy_instance()) {
      ViewT const& image = imagebase.impl();

      // Find pyramid_level and transform image
      TransformRef tx(ResampleTransform(1,1)); // Transform ref doesn't
      int pyramid_level;                       // support a generic construct.
      ImageViewRef<typename ViewT::pixel_type> trans_view = image;
      this->transform_image( input_georef, trans_view, tx, pyramid_level );
      trans_view = per_pixel_filter( trans_view, ClampAlpha() );

      // Calculated affected tiles and print debug statistics
      std::list<TileInfo> tiles;
      this->affected_tiles( bounding_box(image), tx,
                            m_platefile->default_tile_size(),
                            pyramid_level, tiles );

      BBox2i affected_bbox;
      BOOST_FOREACH( TileInfo const& tile, tiles ) {
        affected_bbox.grow( Vector2i(tile.i,tile.j) );
      }
      size_t tiles_size = tiles.size();
      m_platefile->transaction_begin( description, transaction_id_override );

      VW_OUT(InfoMessage, "platefile")
        << "\t    Rasterizing " << tiles_size << " image tiles.\n"
        << "\t    Platefile ID: " << (m_platefile->index_header().platefile_id()) << "\n"
        << "\t    Transaction ID: " << m_platefile->transaction_id() << "\n"
        << "\t    Affected tiles @ root: " << affected_bbox << "\n";

      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      m_platefile->write_request();

      // Add each tile
      progress.report_progress(0);
      BOOST_FOREACH( TileInfo const& tile, tiles ) {
        typedef WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> > Job;

        boost::scoped_ptr<Task> task(
          new Job(m_platefile,
                  tile, pyramid_level,
                  trans_view, false, boost::numeric_cast<int>(tiles_size), progress));
        (*task)();
      }
      progress.report_finished();

      // Sync the index
      m_platefile->sync();

      // Mipmap the tiles.
      if (pyramid_level > 0) {
        std::ostringstream mipmap_str;
        mipmap_str << "\t--> Mipmapping from level " << pyramid_level << ": ";
        this->mipmap(pyramid_level, BBox2i(affected_bbox.min(), affected_bbox.max()+Vector2i(1,1)), m_platefile->transaction_id(),
                     (!tweak_settings_for_terrain), // mipmap preblur = !tweak_settings_for_terrain
                     TerminalProgressCallback( "plate", mipmap_str.str()));
      }

      // Release the blob id lock.
      m_platefile->write_complete();

      // Notify the index that this transaction is complete.  Do not
      // update the read cursor (false).
      m_platefile->transaction_end(false);
    }

  };

  // -------------------------------------------------------------------------
  //                            WRITE PLATEFILE TASK
  // -------------------------------------------------------------------------

  template <class ViewT>
  class WritePlateFileTask : public Task {
    boost::shared_ptr<PlateFile> m_platefile;
    TileInfo m_tile_info;
    int m_level;
    ViewT const& m_view;
    bool m_verbose;
    SubProgressCallback m_progress;

  public:
    WritePlateFileTask(boost::shared_ptr<PlateFile> platefile,
                       TileInfo const& tile_info,
                       int level, ImageViewBase<ViewT> const& view,
                       bool verbose, int total_num_blocks,
                       const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) : m_platefile(platefile),
      m_tile_info(tile_info), m_level(level), m_view(view.impl()),
      m_verbose(verbose), m_progress(progress_callback,0.0,1.0/float(total_num_blocks)) {}

    virtual ~WritePlateFileTask() {}
    virtual void operator() () {
      VW_OUT(DebugMessage, "platefile") << "\t    Generating tile: [ "
                                        << m_tile_info.i << " " << m_tile_info.j
                                        << " @ level " <<  m_level << "]    BBox: "
                                        << m_tile_info.bbox << "\n";

      // Generate the tile from the image data
      ImageView<typename ViewT::pixel_type> tile = crop(m_view, m_tile_info.bbox);

      // If this tile contains no data at all, then we bail early without
      // doing anything.
      if (is_transparent(tile)) {
        m_progress.report_incremental_progress(1.0);
        return;
      }

      // If the data is valid, we write it to the platefile.
      //
      // TODO: This is where we could strip the tile of its alpha
      // channel to save space in the placefile.  This will require a
      // view that strips off the alpha channel.
      //      m_platefile->write_request();

      switch(m_platefile->pixel_format()) {
      case VW_PIXEL_GRAYA:
        switch(m_platefile->channel_type()) {
        case VW_CHANNEL_UINT8:
        case VW_CHANNEL_UINT16:
          m_platefile->write_update(pixel_cast<PixelGrayA<uint8> >(tile), m_tile_info.i, m_tile_info.j, m_level);
          break;
        case VW_CHANNEL_INT16:
          m_platefile->write_update(pixel_cast<PixelGrayA<int16> >(tile), m_tile_info.i, m_tile_info.j, m_level);
          break;
        case VW_CHANNEL_FLOAT32:
          m_platefile->write_update(pixel_cast<PixelGrayA<float32> >(tile), m_tile_info.i, m_tile_info.j, m_level);
          break;
        default:
          vw_throw(NoImplErr() << "Unsupported GrayA channel type in PlateManager.");
        }
        break;
      case VW_PIXEL_RGBA:
        switch(m_platefile->channel_type()) {
        case VW_CHANNEL_UINT8:
          m_platefile->write_update(pixel_cast<PixelRGBA<uint8> >(tile), m_tile_info.i, m_tile_info.j, m_level);
          break;
        default:
          vw_throw(NoImplErr() << "Unsupported RGBA channel type in PlateManager.");
        }
        break;
      default:
        vw_throw(NoImplErr() << "Unsupported pixel type in PlateManager.");
      }

      //      m_platefile->write_complete();
      m_progress.report_incremental_progress(1.0);
    }
  };



}} // namespace vw::plate

#endif // __VW_PLATE_PLATEMANAGER_H__
