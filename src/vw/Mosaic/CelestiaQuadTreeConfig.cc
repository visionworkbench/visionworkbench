// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
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
