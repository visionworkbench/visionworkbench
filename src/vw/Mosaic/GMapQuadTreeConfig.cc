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

  boost::shared_ptr<ImageResource> GMapQuadTreeConfig::tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    return boost::shared_ptr<ImageResource>( DiskImageResource::create( info.filepath, format, info.filetype ) );
  }

  void GMapQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_cull_images( true );
    qtree.set_file_type( "auto" );
    qtree.set_image_path_func( &image_path );
    qtree.set_tile_resource_func( &tile_resource );
  }

} // namespace mosaic
} // namespace vw
