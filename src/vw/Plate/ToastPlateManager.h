// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ImageComposite.h>
#include <vw/Plate/PlateManager.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

// Protocol Buffer
#include <vw/Plate/ProtoBuffers.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------
  
  class ToastPlateManager : public PlateManager {

  public:
  
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      PlateManager(platefile, num_threads) {}

    /// Destructor
    virtual ~ToastPlateManager() {}

    /// Compute the bounding boxes for [ tile_size x tile_size ] tiles
    /// to generate for an image with image_bbox in a TOAST
    /// projection space that is [ resolution x resolution ] pixels in
    /// size.
    ///
    /// NOTE: Tiles in the TOAST projection overlap with their
    /// neighbors ( the last row and last column of pixels is the same
    /// as the top row and left column of the next images.)  This
    /// means that these bounding boxes are a little funny.  This
    /// function computes those bounding boxes for the tiles at the
    /// bottom of the pyramid so that there is the proper amount of
    /// overlap.
    ///
    /// TODO: This function is slow and not very smart -- it computes
    /// all possible bounding boxes before selecting down to the
    /// possible handful that overlap with the image_bbox.  It could
    /// be made to go MUCH faster by being smart about its search.
    std::vector<TileInfo> wwt_image_tiles( BBox2i const& image_bbox, 
                                              int32 const resolution,
                                              int32 const tile_size);

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, cartography::GeoReference const& georef) {

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
      double delta = sqrt(pow(p1.y()-p0.y(),2)+pow(p2.y()-p0.y(),2));
      int pyramid_level = (int)round(log(360/delta/(tile_size-1)) /
                                     log(2));
      
      // Compute the resolution of the TOAST output space at the given
      // pyramid_level.  The formula below was carefully chosen to
      // ensure that the output space at each level is sized just
      // slightly shy of an even power of two; therby allowing for the
      // slight overlap of one line of pixels on the right and the
      // bottom of a toast tile.  This overlap is necessary for proper
      // positioning when these tiles are rendered as texture in a 3D
      // graphics environment.
      int32 resolution = (1<<pyramid_level)*(tile_size-1)+1;

      // Set up the toast transform and compute the bounding box of this
      // image in the toast projection space.
      cartography::ToastTransform toast_tx( georef, resolution );
      BBox2i input_bbox = BBox2i(0,0,image.impl().cols(),image.impl().rows());
      BBox2i output_bbox = toast_tx.forward_bbox(input_bbox);

      std::cout << "\t--> Placing image at level " << pyramid_level 
                << " with bbox " << output_bbox << "\n"
                << "\t    (Total TOAST resolution at this level =  " 
                << resolution << " pixels.)\n";

      // Create the output view and crop it to the proper size.
      ImageViewRef<typename ViewT::pixel_type> toast_view = 
        transform(image,toast_tx, ZeroEdgeExtension(),BicubicInterpolation());

      if( georef.proj4_str()=="+proj=longlat" &&
          fabs(georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(180,0)).x() - image.impl().cols()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.impl().rows()) < 1 ) {
        vw_out(0) << "\t--> Detected global overlay.  " 
                  << "Using cylindrical edge extension to hide the seam.\n";
        toast_view = transform(image,toast_tx,
                               CylindricalEdgeExtension(),BicubicInterpolation());
      } 

      // chop up the image into small chunks
      std::vector<TileInfo> tiles = wwt_image_tiles( output_bbox, resolution,
                                                     m_platefile->default_tile_size());
      
      // And save each tile to the PlateFile
      std::cout << "\t--> Rasterizing " << tiles.size() << " image tiles.\n";
      for (int i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            toast_view)));
      }
      m_queue.join_all();
    }


    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    virtual ImageView<PixelRGBA<uint8> > load_tile( int32 level, int32 x, int32 y );

  };


}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
