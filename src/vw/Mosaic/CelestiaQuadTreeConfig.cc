// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Mosaic/CelestiaQuadTreeConfig.h>

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace mosaic {

  std::string CelestiaQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name(), fs::native );

    Vector2i pos(0,0);
    for ( int i=0; i<(int)name.length(); ++i ) {
      pos *= 2;

      if( name[i]=='2' )      pos += Vector2i(0,1);
      else if( name[i]=='3' ) pos += Vector2i(1,1);
      else if( name[i]=='1' ) pos += Vector2i(1,0);
      else if( name[i]!='0' ) {
        vw_throw(LogicErr() << "Celestia output format incompatible with non-standard quadtree structure");
      }
    }

    int max_val = int(::pow(2, name.length())) >> 1;

    std::ostringstream oss;
    if (name.length() == 0) {
      oss << "original";
    } else {
      oss << "level" << name.length()-1 << "/" << "tx_" << pos.x() << "_" << pos.y()-max_val;
    }

    path /= oss.str();

    return path.native_file_string();
  }

  void CelestiaQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_image_path_func( &image_path );
    qtree.set_cull_images( true );
  }

} // namespace mosaic
} // namespace vw
