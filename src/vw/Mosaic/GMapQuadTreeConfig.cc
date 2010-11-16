// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Mosaic/GMapQuadTreeConfig.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/foreach.hpp>
namespace fs = boost::filesystem;

#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace mosaic {

  // for gmaps, origin is upper left, advancing right and down.
  std::string GMapQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name(), fs::native );

    Vector<size_t, 2> pos(0,0);
    BOOST_FOREACH(char n, name) {
      pos *= 2;
      // this is a pretty weird layout, but it's the one that matches the
      // current nasamaps.js. Don't blame me, I didn't vote for him.
      switch (n) {
        case '0': pos += Vector2i(0,1); break; // Upper Left
        case '1': pos += Vector2i(1,1); break; // Upper Right
        case '2': pos += Vector2i(0,0); break; // Lower Left
        case '3': pos += Vector2i(1,0); break; // Lower Right
        default:
          vw_throw( LogicErr() << "GMap output format incompatible with non-standard quadtree structure" );
      }
    }
    std::ostringstream oss;
    oss << name.length() << "/" << pos.x() << "/" << pos.y();
    path /= oss.str();

    return path.native_file_string();
  }

  boost::shared_ptr<DstImageResource> GMapQuadTreeConfig::tile_resource( QuadTreeGenerator const& /*qtree*/, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    return boost::shared_ptr<DstImageResource>( DiskImageResource::create( info.filepath + info.filetype, format, info.filetype ) );
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

    const double RADIUS = 360. / (2. * M_PI);

    // GMercatorProjection is Mercator. It uses a spherical datum. We choose an
    // arbitrary radius.  This arbitrary radius happens to map mercator onto
    // (roughly) degrees, because image2qtree assumes that.
    cartography::Datum d("Google Maps", "EARTH", "Greenwich", RADIUS, RADIUS, 0);

    cartography::GeoReference r(d);
    r.set_pixel_interpretation(cartography::GeoReference::PixelAsArea);
    r.set_mercator(0,0);

    const double half_circum = M_PI * RADIUS;
    Matrix3x3 transform;
    transform(0,0) = 2.0 * (half_circum / xresolution);
    transform(0,2) = -half_circum;
    transform(1,1) = -2.0 * (half_circum / yresolution);
    transform(1,2) = half_circum;
    transform(2,2) = 1;
    r.set_transform(transform);

    return r;
  }

} // namespace mosaic
} // namespace vw
