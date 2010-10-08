// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Mosaic/GMapQuadTreeConfig.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace mosaic {

  // for gmaps, origin is upper left, advancing right and down.
  std::string GMapQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name(), fs::native );

    Vector2i pos(0,0);
    for ( int i=0; i<(int)name.length(); ++i ) {
      pos *= 2;
      if( name[i]=='1' ) pos += Vector2i(1,0);
      else if( name[i]=='2' ) pos += Vector2i(0,1);
      else if( name[i]=='3' ) pos += Vector2i(1,1);
      else if( name[i]!='0' ) {
        vw_throw( LogicErr() << "GMap output format incompatible with non-standard quadtree structure" );
      }
    }
    std::ostringstream oss;
    oss << name.length() << "/" << pos.x() << "/" << pos.y();
    path /= oss.str();

    return path.native_file_string();
  }

  boost::shared_ptr<ImageResource> GMapQuadTreeConfig::tile_resource( QuadTreeGenerator const& /*qtree*/, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    return boost::shared_ptr<ImageResource>( DiskImageResource::create( info.filepath + info.filetype, format, info.filetype ) );
  }

  void GMapQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_cull_images( true );
    qtree.set_file_type( "jpg" );
    qtree.set_image_path_func( &image_path );
    qtree.set_tile_resource_func( &tile_resource );
  }

  cartography::GeoReference GMapQuadTreeConfig::output_georef(uint32 xresolution, uint32 yresolution) {
    if (yresolution == 0)
      yresolution = xresolution;

    VW_ASSERT(xresolution == yresolution, LogicErr() << "GMap requires square pixels");

    // GMercatorProjection is Mercator, -180 -> 180
    cartography::GeoReference r;
    r.set_pixel_interpretation(cartography::GeoReference::PixelAsArea);
    r.set_mercator(0,0);

    Matrix3x3 transform;
    transform(0,0) = 360.0 / xresolution;
    transform(0,2) = -180;
    transform(1,1) = -360.0 / yresolution;
    transform(1,2) = 180;
    transform(2,2) = 1;
    r.set_transform(transform);

    return r;
  }

} // namespace mosaic
} // namespace vw
