// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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


#include <vw/Mosaic/QuadTreeGenerator.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace mosaic {

  std::string QuadTreeGenerator::simple_image_path::operator()( QuadTreeGenerator const& qtree, std::string const& name ) {
    fs::path path( qtree.get_name() );
    path /= "r" + name;
    return path.string();
  }

  std::string QuadTreeGenerator::tiered_image_path::operator()( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory ) {
    fs::path path( qtree.get_name() );

    std::string rname = "r" + name;

    for (int32 i= 0; i< (int32)rname.length() - levels_per_directory; i += levels_per_directory) {
      path /= rname.substr( i, levels_per_directory );
    }
    path /= rname;

    return path.string();
  }

  std::string QuadTreeGenerator::named_tiered_image_path::operator()( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory ) {
    fs::path path( qtree.get_name() );

    if( name.length() == 0 ) {
      path /= fs::path(path).replace_extension("").filename();
    }
    else {
      for ( int32 i=0; i<(int32)name.length() - levels_per_directory; i+=levels_per_directory ) {
        path /= name.substr( i, levels_per_directory );
      }
      path /= name;
    }

    return path.string();
  }

  std::vector<std::pair<std::string,vw::BBox2i> > QuadTreeGenerator::default_branch_func::operator()( QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region ) {
    std::vector<std::pair<std::string,vw::BBox2i> > children;
    if( region.width() > qtree.get_tile_size() && region.height() > qtree.get_tile_size() ) {
      children.push_back( std::make_pair( name + "0", BBox2i( (region + region.min()) / 2 ) ) );
      children.push_back( std::make_pair( name + "1", BBox2i( (region + Vector2i(region.max().x(),region.min().y())) / 2 ) ) );
      children.push_back( std::make_pair( name + "2", BBox2i( (region + Vector2i(region.min().x(),region.max().y())) / 2 ) ) );
      children.push_back( std::make_pair( name + "3", BBox2i( (region + region.max()) / 2 ) ) );
    }
    return children;
  }

  boost::shared_ptr<DstImageResource> QuadTreeGenerator::default_tile_resource_func::operator()( QuadTreeGenerator const&, TileInfo const& info, ImageFormat const& format ) {
    create_directories( fs::path( info.filepath ).parent_path() );
    return boost::shared_ptr<DstImageResource>( DiskImageResource::create( info.filepath+info.filetype, format ) );
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
