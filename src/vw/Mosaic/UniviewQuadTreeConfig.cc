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


#include <vw/Mosaic/UniviewQuadTreeConfig.h>

#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <functional>
namespace fs = boost::filesystem;

#include <vw/FileIO/DiskImageResourcePNG.h>

namespace vw {
namespace mosaic {

  // A special subclass of DiskImageResourcePNG that produces the bizarro
  // 16-bit PNG files that Uniview uses for terrain data.  They're standard
  // uint16 image files, except the data is interpreted as signed int16
  // data instead.
  class UniviewTerrainResource : public DiskImageResourcePNG {

    // A convenience function to force the format hint that we
    // pass to DiskImageResourcePNG to be single-channel uint16
    static ImageFormat make_uint16( ImageFormat const& format ) {
      ImageFormat result = format;
      result.pixel_format = VW_PIXEL_GRAY;
      result.channel_type = VW_CHANNEL_UINT16;
      return result;
    }

  public:
    // We only provide a constructor for image *creation*, at least for now.
    UniviewTerrainResource( std::string const& filename, ImageFormat const& format )
      : DiskImageResourcePNG( filename, make_uint16(format) )
    {}

    // First we convert to single-channel signed int16, then we spoof that as
    // uint16 data and pass it along to DiskImageResourcePNG to write.
    void write( ImageBuffer const& src, BBox2i const& bbox ) {
      ImageView<PixelGray<int16> > im_buf( src.format.cols, src.format.rows );
      ImageBuffer buffer = im_buf.buffer();
      convert( buffer, src );
      buffer.format.channel_type = VW_CHANNEL_UINT16;
      DiskImageResourcePNG::write( buffer, bbox );
    }
  };

  void UniviewQuadTreeConfig::set_module(const std::string& module) {
    m_module_name = module;
  }

  std::string UniviewQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name() );

    Vector2i pos(0,0);
    for ( int i=0; i<(int)name.length(); ++i ) {
      pos *= 2;
      if( name[i]=='0' ) pos += Vector2i(0,1);
      else if( name[i]=='1' ) pos += Vector2i(1,1);
      else if( name[i]=='3' ) pos += Vector2i(1,0);
      else if( name[i]!='2' ) {
        vw_throw(LogicErr() << "Uniview output format incompatible with non-standard quadtree structure");
      }
    }
    std::ostringstream oss;
    if (name.length() == 0)
      oss << "global";
    else
      oss << name.length()-1 << "/" << pos.y() << "/" << pos.x();
    path /= oss.str();

    return path.string();
  }


  boost::shared_ptr<DstImageResource> UniviewQuadTreeConfig::terrain_tile_resource( QuadTreeGenerator const& /*qtree*/,QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) {
    fs::create_directories(fs::path( info.filepath ).parent_path());
    return boost::shared_ptr<DstImageResource>( new UniviewTerrainResource( info.filepath+info.filetype, format ) );
  }


  void UniviewQuadTreeConfig::configure( QuadTreeGenerator &qtree ) const {
    qtree.set_image_path_func( &image_path );
    if( m_terrain ) qtree.set_tile_resource_func( &terrain_tile_resource );
    qtree.set_metadata_func( boost::bind(&UniviewQuadTreeConfig::metadata_func,this,
                                         std::placeholders::_1, std::placeholders::_2));
  }

  // TODO: Is this actually the right function for Uniview?
  cartography::GeoReference UniviewQuadTreeConfig::output_georef(uint32 xresolution, uint32 yresolution) {
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

  void UniviewQuadTreeConfig::metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const {
    // root node
    if (info.name.empty()) {
      fs::path out_path(qtree.get_name());
      fs::path conf_path = fs::path(out_path).replace_extension(".conf");

      fs::ofstream conf( conf_path );
      if(m_terrain) {
        conf << "// Terrain\n";
        conf << "HeightmapCacheLocation=modules/" << m_module_name << "/Offlinedatasets/" << out_path.filename() << "/Terrain/\n";
        conf << "HeightmapCallstring=Generated by the NASA Vision Workbench image2qtree tool.\n";
        conf << "HeightmapFormat=" << qtree.get_file_type() << '\n';
        conf << "NrHeightmapLevels=" << qtree.get_tree_levels()-1 << '\n';
        conf << "NrLevelsPerHeightmap=1\n";
      } else {
        conf << "[Offlinedataset]\n";
        conf << "NrRows=1\n";
        conf << "NrColumns=2\n";
        conf << "Bbox= -180 -90 180 90\n";
        conf << "DatasetTitle=" << out_path.filename() << "\n";
        conf << "Tessellation=19\n\n";

        conf << "// Texture\n";
        conf << "TextureCacheLocation=modules/" << m_module_name << "/Offlinedatasets/" << out_path.filename() << "/Texture/\n";
        conf << "TextureCallstring=Generated by the NASA Vision Workbench image2qtree tool.\n";
        conf << "TextureFormat=" << qtree.get_file_type() << "\n";
        conf << "TextureLevels= " << qtree.get_tree_levels()-1 << "\n";
        conf << "TextureSize= " << qtree.get_tile_size() << "\n\n";
      }
      std::cout << "\nNote: You must merge the texture and terrain config files into a single file (Terrain info should go below texture info.)" << std::endl;
      std::cout << "Both output sets should be in the same directory, with the texture in a subdirectory named Texture and the terrain in a subdirectory named Terrain." << std::endl;
    }
  }

} // namespace mosaic
} // namespace vw
