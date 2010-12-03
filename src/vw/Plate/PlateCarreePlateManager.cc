// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Image/Filter.h>

using namespace vw::platefile;
using namespace vw;

template <class PixelT>
void PlateCarreePlateManager<PixelT>::affected_tiles(
                          BBox2i const& image_size,
                          TransformRef const& tx, int tile_size,
                          int level, std::list<TileInfo>& tiles ) const {
  PlateManager<PixelT>::affected_tiles( image_size, tx, tile_size,
                                        level, tiles );

  int32 level_size_tiles = 1 << level;
  int32 level_size_pixels = level_size_tiles*tile_size;

  // Correct tiles that cross boundaries
  BOOST_FOREACH( TileInfo& tile, tiles ) {
    if (tile.i < 0) {
      tile.i += level_size_tiles;
      tile.bbox.min() += Vector2i(level_size_pixels,0);
      tile.bbox.max() += Vector2i(level_size_pixels,0);
    }
    if (tile.i >= level_size_tiles) {
      tile.i -= level_size_tiles;
      tile.bbox.min() -= Vector2i(level_size_pixels,0);
      tile.bbox.max() -= Vector2i(level_size_pixels,0);
    }
  }
}

template <class PixelT>
void PlateCarreePlateManager<PixelT>::transform_image(
                          cartography::GeoReference const& georef,
                          ImageViewRef<PixelT>& image,
                          TransformRef& txref, int& level ) const {
  // First .. correct input transform if need be
  cartography::GeoReference input_georef = georef;
  if ( boost::contains(input_georef.proj4_str(),"+proj=longlat") ) {
    Matrix3x3 transform = input_georef.transform();
    // Correct if it is so far right it is not visible.
    if ( transform(0,2) > 180 )
      transform(0,2) -= 360;
    if ( input_georef.pixel_to_lonlat(Vector2(image.cols()-1,0))[0] < -180.0 )
      transform(0,2) += 360;
    input_georef.set_transform(transform);
  }

  // Create temporary transform to work out the resolution
  cartography::GeoReference output_georef;
  output_georef.set_datum(input_georef.datum());
  int resolution = 256;
  cartography::GeoTransform geotx( input_georef, output_georef );
  // Calculate the best resolution at 5 different points in the image
  Vector2 res_pixel[5];
  res_pixel[0] = Vector2( image.cols()/2, image.rows()/2 );
  res_pixel[1] = Vector2( image.cols()/2 + image.cols()/4,
                          image.rows()/2 );
  res_pixel[2] = Vector2( image.cols()/2 - image.cols()/4,
                          image.rows()/2 );
  res_pixel[3] = Vector2( image.cols()/2,
                          image.rows()/2 + image.rows()/4 );
  res_pixel[4] = Vector2( image.cols()/2,
                          image.rows()/2 - image.rows()/4 );
  int res;
  for(int i=0; i < 5; i++) {
    res = cartography::output::kml::compute_resolution(geotx, res_pixel[i]);
    if( res > resolution ) resolution = res;
  }

  // Round the resolution to the nearest power of two.  The
  // base of the pyramid is 2^8 or 256x256 pixels.
  level = static_cast<int>(ceil(log(resolution) / log(2))) - 8;
  output_georef = this->georeference( level );

  // Rebuild the transform
  geotx = cartography::GeoTransform( input_georef, output_georef );
  BBox2i output_bbox = geotx.forward_bbox(bounding_box(image));
  vw_out() << "\t    Placing image at level " << level
           << " with bbox " << output_bbox << "\n"
           << "\t    (Total KML resolution at this level =  "
           << resolution << " pixels.)\n";

  // Perform transform and rewrite to input
  ImageViewRef<PixelT> holding =
    transform( image, geotx, ZeroEdgeExtension(),
               BicubicInterpolation() );
  if ( output_bbox.max().x() > resolution ) {
    // Determine if we are too far east
    image = edge_extend(crop(holding,resolution/2,0,resolution,resolution),
                        -resolution/2,0,resolution,resolution,
                        PeriodicEdgeExtension());
  } else if ( output_bbox.min().x() < 0 ) {
    // Determine if we are too far west
    image = edge_extend(crop(holding,-resolution/2,0,resolution,resolution),
                        resolution/2,0,resolution,resolution,
                        PeriodicEdgeExtension());
  } else {
    image = holding;
  }
  txref = TransformRef(geotx);
}

template <class PixelT>
cartography::GeoReference
PlateCarreePlateManager<PixelT>::georeference( int level ) const {
  int tile_size = this->m_platefile->default_tile_size();
  int resolution = (1<<level)*tile_size;

  cartography::GeoReference r;
  r.set_pixel_interpretation(cartography::GeoReference::PixelAsArea);

  // Set projecion space to be between -180 and 180.
  Matrix3x3 transform;
  transform(0,0) = 360.0 / resolution;
  transform(0,2) = -180.0;
  transform(1,1) = -360.0 / resolution;
  transform(1,2) = 180.0;
  transform(2,2) = 1;
  r.set_transform(transform);

  return r;
}

// Explicit template instantiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_PLATE_CARREE_TYPES(PIXELT)                            \
  template void                                                              \
  PlateCarreePlateManager<PIXELT >::transform_image(                         \
                                    cartography::GeoReference const& georef, \
                                    ImageViewRef<PIXELT >& image,            \
                                    TransformRef& txref, int& level ) const; \
  template void                                                              \
  PlateCarreePlateManager<PIXELT >::affected_tiles(                          \
                                    BBox2i const& image_size,                \
                                    TransformRef const& tx, int tile_size,   \
                              int level, std::list<TileInfo>& tiles ) const; \
  template cartography::GeoReference                                         \
  PlateCarreePlateManager<PIXELT >::georeference( int level ) const;

  VW_INSTANTIATE_PLATE_CARREE_TYPES(PixelGrayA<uint8>)
  VW_INSTANTIATE_PLATE_CARREE_TYPES(PixelGrayA<int16>)
  VW_INSTANTIATE_PLATE_CARREE_TYPES(PixelGrayA<float32>)
  VW_INSTANTIATE_PLATE_CARREE_TYPES(PixelRGBA<uint8>)

}}
