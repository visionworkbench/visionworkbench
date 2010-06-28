// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/PlateManager.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <vector>
#include <fstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace platefile {
    // Compute the bounding boxes for [ tile_size x tile_size ] tiles
    // to generate for an image with image_bbox in a TOAST
    // projection space that is [ resolution x resolution ] pixels in
    // size.

    // NOTE: Tiles in the TOAST projection overlap with their
    // neighbors ( the last row and last column of pixels is the same
    // as the top row and left column of the next images.)  This
    // means that these bounding boxes are a little funny.  This
    // function computes those bounding boxes for the tiles at the
    // bottom of the pyramid so that there is the proper amount of
    // overlap.

    // TODO: This function is slow and not very smart -- it computes
    // all possible bounding boxes before selecting down to the
    // possible handful that overlap with the tile_bbox.  It could
    // be made to go MUCH faster by being smart about its search.
    std::vector<TileInfo> wwt_image_tiles( BBox2i const& input_bbox,
                                           cartography::ToastTransform const& toast_tx,
                                           BBox2i const& tile_bbox,
                                           int32 const resolution,
                                           int32 const tile_size);


  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------

  template <class PixelT>
  class ToastPlateManager : public PlateManager {

    int32 m_resolution;

    // Private variables
    FifoWorkQueue m_queue;

    struct CacheEntry {
      int32 level, x, y, transaction_id;
      ImageView<PixelT> tile;
    };

    // The cache is declared mutable because it is modified by the
    // otherwise const fetch_child_tile() method.
    typedef std::list<CacheEntry> cache_t;
    mutable cache_t m_cache;  

 public:
  
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile) : 
      PlateManager(platefile), m_queue(1) {} // Use 1 thread for now...

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, std::string const& description, 
                int transaction_id_override, cartography::GeoReference const& georef, 
                bool tweak_settings_for_terrain, bool verbose = false,
                const ProgressCallback &progress = ProgressCallback::dummy_instance()) {
      
      // Compute the pyramid level at which to store this image.  The
      // number of required levels is broken down so that the very top
      // of the pyramid covers the entire globe and has a size of
      // tile_size.
      int tile_size = m_platefile->default_tile_size();
      Vector2 p0 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2,
                                                  image.impl().rows()/2));
      Vector2 p1 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2+1,
                                                  image.impl().rows()/2));
      Vector2 p2 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2,
                                                  image.impl().rows()/2+1));
      double degrees_per_pixel = sqrt(pow(p1.y()-p0.y(),2)+pow(p2.y()-p0.y(),2));
      int pyramid_level = (int)round(log(360/degrees_per_pixel/(tile_size-1)) / log(2));
      
      // Compute the resolution of the TOAST output space at the given
      // pyramid_level.  The formula below was carefully chosen to
      // ensure that the output space at each level is sized just
      // slightly shy of an even power of two; therby allowing for the
      // slight overlap of one line of pixels on the right and the
      // bottom of a toast tile.  This overlap is necessary for proper
      // positioning when these tiles are rendered as texture in a 3D
      // graphics environment.
      m_resolution = (1<<pyramid_level)*(tile_size-1)+1;

      // Set up the toast transform and compute the bounding box of this
      // image in the toast projection space.
      
      cartography::ToastTransform toast_tx( georef, m_resolution );
      BBox2i input_bbox = BBox2i(0,0,image.impl().cols(),image.impl().rows());
      BBox2i output_bbox = toast_tx.forward_bbox(input_bbox);

      vw_out(InfoMessage, "platefile")
        << "\t    Placing image at level " << pyramid_level
        << " with bbox " << output_bbox << "\n"
        << "\t    (Total TOAST resolution at this level =  "
        << m_resolution << " pixels.)\n";

      // Create the output view and crop it to the proper size.
      ImageViewRef<typename ViewT::pixel_type> toast_view;
      // if (tweak_settings_for_terrain) {
      //     toast_view = transform(image,toast_tx, ZeroEdgeExtension(), BilinearInterpolation());
      // } else {
          toast_view = transform(image,toast_tx, ZeroEdgeExtension(), BicubicInterpolation());
      // }
      if( (boost::trim_copy(georef.proj4_str())=="+proj=longlat") &&
          (fabs(georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1) &&
          (fabs(georef.lonlat_to_pixel(Vector2(180,0)).x() - image.impl().cols()) < 1) &&
          (fabs(georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1) &&
          (fabs(georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.impl().rows()) < 1) ) {
        vw_out() << "\t--> Detected global overlay.  " 
                  << "Using cylindrical edge extension to hide the seam.\n";
        // if (tweak_settings_for_terrain) {
        //   toast_view = transform(image,toast_tx,
        //                          CylindricalEdgeExtension(),BilinearInterpolation());
        // } else {
          toast_view = transform(image,toast_tx,
                                 CylindricalEdgeExtension(),BicubicInterpolation());
          //        }
      } 

      // chop up the image into small chunks
      std::vector<TileInfo> tiles = wwt_image_tiles( input_bbox, toast_tx, output_bbox,
                                                     m_resolution,
                                                     m_platefile->default_tile_size());

      // Compute the affected tiles.
      BBox2i affected_tiles_bbox;
      for (size_t i = 0; i < tiles.size(); ++i) 
        affected_tiles_bbox.grow(Vector2i(tiles[i].i,tiles[i].j));

      // Obtain a transaction ID for this image.  
      //
      // Note: the user may have specified a transaction_id to use,
      // which was passed in with transaction_id_override.  If not,
      // then transaction_id_override == -1, and we get an
      // automatically assigned t_id.
      int transaction_id = m_platefile->transaction_request(description, 
                                                            transaction_id_override);

      vw_out(InfoMessage, "platefile")
        << "\t    Rasterizing " << tiles.size() << " image tiles.\n"
        << "\t    Platefile ID: " << (m_platefile->index_header().platefile_id()) << "\n"
        << "\t    Transaction ID: " << transaction_id << "\n"
        << "\t    Affected tiles @ root: " << affected_tiles_bbox << "\n";

      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      m_platefile->write_request();

      // Add each tile.
      progress.report_progress(0);
      for (size_t i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            transaction_id,
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            toast_view, 
                                                                            tweak_settings_for_terrain,
                                                                            verbose, // verbose
                                                                            tiles.size(),
                                                                            progress)));
      }
      m_queue.join_all();
      progress.report_finished();
      
      // Sync the index
      m_platefile->sync();

      // Mipmap the tiles.
      if (m_platefile->num_levels() > 1) {
        std::ostringstream mipmap_str;
        mipmap_str << "\t--> Mipmapping from level " << pyramid_level << ": ";
        this->mipmap(pyramid_level, affected_tiles_bbox, transaction_id, 
                     !tweak_settings_for_terrain, // mipmap preblur = !tweak_settings_for_terrain
                     TerminalProgressCallback( "plate", mipmap_str.str()));
      }

      // Release the blob id lock.
      m_platefile->write_complete();

      // Notify the index that this transaction is complete.  Do not
      // update the read cursor (false).
      m_platefile->transaction_complete(transaction_id, false);
    }

    /// This function generates a specific mipmap tile at the given
    /// col, row, and level, and transaction_id.  
    virtual void generate_mipmap_tile(int col, int row, int level, int transaction_id, bool preblur) const;

    ImageView<PixelT> fetch_child_tile(int x, int y, int level, int transaction_id) const;
  };


}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
