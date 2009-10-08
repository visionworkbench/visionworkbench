// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_KML_PLATEMANAGER_H__
#define __VW_PLATE_KML_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
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
  
  class KmlPlateManager : public PlateManager {

  public:
  
    KmlPlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      PlateManager(platefile, num_threads) {}

    /// Destructor
    virtual ~KmlPlateManager() {}

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
    std::vector<TileInfo> kml_image_tiles( BBox2i const& image_bbox, 
                                           int32 const resolution,
                                           int32 const tile_size);

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, 
                cartography::GeoReference const& input_georef) {

      
      // We build the ouput georeference from the input image's
      // georeference by copying the datum.
      cartography::GeoReference output_georef;
      output_georef.set_datum(input_georef.datum());
      int resolution = 1024;

      // Right now, we only need a WGS84 output geoereference to compute 
      // the resolution. The rest of the output info will get set later.
      cartography::GeoTransform geotx( input_georef, output_georef );
     
      // Calculate the best resolution at 5 different points in the image, 
      // as occasionally there's a singularity at the center pixel that 
      // makes it extremely tiny (such as in pole-centered images).
      const int cols = image.impl().cols();
      const int rows = image.impl().rows();
      Vector2 res_pixel[5];
      res_pixel[0] = Vector2( cols/2, rows/2 );
      res_pixel[1] = Vector2( cols/2 + cols/4, rows/2 );
      res_pixel[2] = Vector2( cols/2 - cols/4, rows/2 );
      res_pixel[3] = Vector2( cols/2, rows/2 + rows/4 );
      res_pixel[4] = Vector2 (cols/2, rows/2 - rows/4 );
      int res;
      for(int i=0; i < 5; i++) {
        res = vw::cartography::output::kml::compute_resolution(geotx, res_pixel[i]);
        if( res > resolution ) resolution = res;
      }
        
      // Round the resolution to the nearest power of two
      int pyramid_level = (int)ceil(log(resolution) / log(2));
      int tile_size = m_platefile->default_tile_size();
      resolution = (1<<pyramid_level)*tile_size;

      // Set up the KML transform and compute the bounding box of this
      // image in the KML projection space.
      cartography::GeoTransform kml_tx( input_georef, output_georef );
      BBox2i input_bbox = BBox2i(0,0,image.impl().cols(),image.impl().rows());
      BBox2i bbox = kml_tx.forward_bbox(input_bbox);
      std::cout << "\t--> Tx bbox: " << bbox << "\n";

      // Compute a tighter Google Earth coordinate system aligned bounding box.
      bbox.crop( BBox2i(0,0,resolution,resolution) );
      BBox2i output_bbox = BBox2i( (bbox.min().x()/resolution)*resolution, 
                                   (bbox.min().y()/resolution)*resolution, 
                                   resolution, resolution );
      if( ! output_bbox.contains( bbox ) ) {
        if( output_bbox.max().x() == resolution ) output_bbox.min().x() -= resolution;
        else output_bbox.max().x() += resolution;
        if( output_bbox.max().y() == resolution ) output_bbox.min().y() -= resolution;
        else output_bbox.max().y() += resolution;
      }
      
      BBox2i ll_bbox = BBox2( -180.0 + (360.0*output_bbox.min().x())/resolution,
                              180.0 - (360.0*output_bbox.max().y())/resolution,
                              (360.0*output_bbox.width())/resolution,
                              (360.0*output_bbox.height())/resolution );
      
      std::cout << "\t--> Placing image at level " << pyramid_level 
                << " with bbox " << output_bbox << "\n"
                << "\t    (Total KML resolution at this level =  " 
                << resolution << " pixels.)\n";
      std::cout << "\t--> Bounding box in geographic space: " << ll_bbox << "\n";
      
      // Create the output view and crop it to the proper size.
      ImageViewRef<typename ViewT::pixel_type> kml_view = 
        transform(image,kml_tx, ZeroEdgeExtension(),BicubicInterpolation());
      
      if( input_georef.proj4_str()=="+proj=longlat" &&
          fabs(input_georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(180,0)).x() - image.impl().cols()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.impl().rows()) < 1 ) {
        vw_out(0) << "\t--> Detected global overlay.  " 
                  << "Using cylindrical edge extension to hide the seam.\n";
        kml_view = transform(image,kml_tx,
                             CylindricalEdgeExtension(), BicubicInterpolation());
      } 

      // chop up the image into small chunks
      std::vector<TileInfo> tiles = kml_image_tiles( output_bbox, resolution,
                                                     m_platefile->default_tile_size());


      // And save each tile to the PlateFile
      std::cout << "\t--> Rasterizing " << tiles.size() << " image tiles.\n";
      for (int i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            kml_view)));
      }
      m_queue.join_all();
    }


    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    virtual ImageView<PixelRGBA<uint8> > load_tile( int32 level, int32 x, int32 y );
  };


}} // namespace vw::plate

#endif // __VW_PLATE_KML_PLATEMANAGER_H__
