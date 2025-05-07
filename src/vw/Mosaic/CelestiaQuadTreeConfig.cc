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


#include <vw/Mosaic/CelestiaQuadTreeConfig.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <functional>

namespace fs = boost::filesystem;

namespace vw {
namespace mosaic {

  void CelestiaQuadTreeConfig::set_module(const std::string& module) {
    m_module_name = module;
  }

  std::string CelestiaQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name() );

    Vector<size_t,2> pos(0,0);
    for ( size_t i=0; i < name.length(); ++i ) {
      pos *= 2;

      if( name[i]=='2' )      pos += Vector2i(0,1);
      else if( name[i]=='3' ) pos += Vector2i(1,1);
      else if( name[i]=='1' ) pos += Vector2i(1,0);
      else if( name[i]!='0' ) {
        vw_throw(LogicErr() << "Celestia output format incompatible with non-standard quadtree structure");
      }
    }

    std::ostringstream oss;
    if (name.length() == 0) {
      oss << "original";
    } else {
      size_t max_val = (1 << (name.length()-1));
      oss << "level" << name.length()-1 << "/" << "tx_" << pos.x() << "_" << max_val - pos.y();
    }

    path /= oss.str();

    return path.string();
  }

  void CelestiaQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_image_path_func( &image_path );
    qtree.set_cull_images( true );
    qtree.set_metadata_func( boost::bind(&CelestiaQuadTreeConfig::metadata_func,this,
                                         std::placeholders::_1, std::placeholders::_2));
  }

  // TODO: Is this actually the right function for Celestia?
  cartography::GeoReference CelestiaQuadTreeConfig::output_georef(uint32 xresolution, uint32 yresolution) {
    if (yresolution == 0)
      yresolution = xresolution;

    VW_ASSERT(xresolution == yresolution, LogicErr() << "TMS requires square pixels");

    cartography::GeoReference r;
    r.set_pixel_interpretation(cartography::GeoReference::PixelAsArea);

    // Note: the global TMS pixel space extends from +270 to -90
    // latitude, so that the lower-left hand corner is tile-
    // aligned, since TMS uses an origin in the lower left.
    Matrix3x3 transform;
    transform(0,0) = 360.0 / xresolution;
    transform(0,2) = -180;
    transform(1,1) = -360.0 / yresolution;
    transform(1,2) = 270;
    transform(2,2) = 1;
    r.set_transform(transform);

    return r;
  }

  void CelestiaQuadTreeConfig::metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const {
    // root node
    if (info.name.empty()) {
      fs::path out_path(qtree.get_name());
      fs::path ctx_path = fs::path(out_path).replace_extension(".ctx");
      fs::path ssc_path = fs::path(out_path).replace_extension(".ssc");

      // TileSize is a heuristic hint to celestia... we make it one step larger
      // at the suggestion of other people who make VTs (it makes it load
      // higher resolutions sooner, I think)
      fs::ofstream ctx( ctx_path );
      ctx << "VirtualTexture\n";
      ctx << "{\n";
      ctx << "        ImageDirectory \"" << out_path.filename() << "\"\n";
      ctx << "        BaseSplit 0\n";
      ctx << "        TileSize " << (qtree.get_tile_size() >> 1) << "\n";
      ctx << "        TileType \"" << qtree.get_file_type() << "\"\n";
      ctx << "}\n";
      ctx.close();

      fs::ofstream ssc( ssc_path );

      ssc << "AltSurface \"" << out_path.filename() << "\" \"" << m_module_name << "\"\n";
      ssc << "{\n";
      ssc << "    Texture \"" << ctx_path.filename() << "\"\n";
      ssc << "}\n";
      ssc.close();

      std::cout << "\nPlace " << ssc_path << " in Celestia's extras dir" << std::endl;
      std::cout << "Place " << ctx_path << " and the output dir ("
                            << out_path << ") in extras/textures/hires" << std::endl;
    }
  }

}} // namespace vw::mosaic
