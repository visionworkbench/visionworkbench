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

#include <vw/Mosaic/QuadTreeGenerator.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace mosaic {

  std::string QuadTreeGenerator::simple_image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name(), fs::native );
    path /= "r" + name;
    return path.native_file_string();
  }

  std::string QuadTreeGenerator::tiered_image_path( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory ) {
    fs::path path( qtree.get_name(), fs::native );
    
    std::string rname = "r" + name;
    
    for (int32 i= 0; i< (int32)rname.length() - levels_per_directory; i += levels_per_directory) {
      path /= rname.substr( i, levels_per_directory );
    }
    path /= rname;
    
    return path.native_file_string();
  }

  std::string QuadTreeGenerator::named_tiered_image_path( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory ) {
    fs::path path( qtree.get_name(), fs::native );
    
    if( name.length() == 0 ) {
      path /= change_extension( path, "" ).leaf();
    }
    else {
      for ( int32 i=0; i<(int32)name.length() - levels_per_directory; i+=levels_per_directory ) {
	path /= name.substr( i, levels_per_directory );
      }
      path /= name;
    }
    
    return path.native_file_string();
  }

  std::vector<std::pair<std::string,vw::BBox2i> > QuadTreeGenerator::default_branch_func( QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region ) {
    std::vector<std::pair<std::string,vw::BBox2i> > children;
    if( region.width() > qtree.get_tile_size() && region.height() > qtree.get_tile_size() ) {
      children.push_back( std::make_pair( name + "0", BBox2i( (region + region.min()) / 2 ) ) );
      children.push_back( std::make_pair( name + "1", BBox2i( (region + Vector2i(region.max().x(),region.min().y())) / 2 ) ) );
      children.push_back( std::make_pair( name + "2", BBox2i( (region + Vector2i(region.min().x(),region.max().y())) / 2 ) ) );
      children.push_back( std::make_pair( name + "3", BBox2i( (region + region.max()) / 2 ) ) );
    }
    return children;
  }

  boost::shared_ptr<ImageResource> QuadTreeGenerator::default_tile_resource_func( QuadTreeGenerator const&, TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    return boost::shared_ptr<ImageResource>( DiskImageResource::create( info.filepath+info.filetype, format ) );
  }

  void QuadTreeGenerator::generate( const ProgressCallback &progress_callback ) {
    ScopedWatch sw("QuadTreeGenerator::generate");
    int32 tree_levels = get_tree_levels();

    vw_out(DebugMessage, "mosaic") << "Using tile size: " << m_tile_size << " pixels" << std::endl;
    vw_out(DebugMessage, "mosaic") << "Generating tile files of type: " << m_file_type << std::endl;
    vw_out(DebugMessage, "mosaic") << "Generating quadtree with " << tree_levels << " levels." << std::endl;

    BBox2i region_bbox = BBox2i(0,0,m_tile_size,m_tile_size) * (1<<(tree_levels-1));
    m_processor->generate( region_bbox, progress_callback );

    progress_callback.report_finished();
  }

} // namespacw vw::mosaic
} // namespace vw
