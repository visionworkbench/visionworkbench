// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Mosaic/UniviewQuadTreeConfig.h>

#include <sstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
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
  

  std::string UniviewQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name(), fs::native );

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
    
    return path.native_file_string();
  }


  boost::shared_ptr<ImageResource> UniviewQuadTreeConfig::terrain_tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    return boost::shared_ptr<ImageResource>( new UniviewTerrainResource( info.filepath+info.filetype, format ) );
  }

  
  void UniviewQuadTreeConfig::configure( QuadTreeGenerator &qtree ) const {
    qtree.set_image_path_func( &image_path );
    if( m_terrain ) qtree.set_tile_resource_func( &terrain_tile_resource );
  }

} // namespace mosaic
} // namespace vw
