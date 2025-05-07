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


#include <vw/Mosaic/GigapanQuadTreeConfig.h>

#include <boost/bind.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace mosaic {

  struct GigapanQuadTreeConfigData {
    BBox2 m_longlat_bbox;

    std::vector<std::pair<std::string,vw::BBox2i> > branch_func( QuadTreeGenerator const&, 
                                                                 std::string       const& name, 
                                                                 BBox2i            const& region ) const;
    void metadata_func( QuadTreeGenerator const&, QuadTreeGenerator::TileInfo const& info ) const;

    public:
      GigapanQuadTreeConfigData() : m_longlat_bbox(0, 0, 0, 0) {}
  }; // End struct

  GigapanQuadTreeConfig::GigapanQuadTreeConfig()
    : m_data( new GigapanQuadTreeConfigData() )
  {}

  std::string GigapanQuadTreeConfig::image_path( QuadTreeGenerator const& qtree, std::string const& name ) {
    return QuadTreeGenerator::tiered_image_path()(qtree, name);
  }

  void GigapanQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_image_path_func( &image_path );
    qtree.set_cull_images( true );
    qtree.set_metadata_func( boost::bind(&GigapanQuadTreeConfigData::metadata_func,m_data,_1,_2) );

    if (m_data->m_longlat_bbox.width() != 0 || m_data->m_longlat_bbox.height() != 0) {
      qtree.set_branch_func( boost::bind(&GigapanQuadTreeConfigData::branch_func,m_data,_1,_2,_3) );
    }
  }

  void GigapanQuadTreeConfig::set_longlat_bbox( BBox2 const& bbox ) {
    m_data->m_longlat_bbox = bbox;
  }

  std::vector<std::pair<std::string,vw::BBox2i> > GigapanQuadTreeConfigData::branch_func( QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region ) const {
    std::vector<std::pair<std::string,vw::BBox2i> > children;
    
    if( region.height() <= qtree.get_tile_size() )
      return children;
      
    Vector2i dims = qtree.get_dimensions();
    double aspect_ratio = 2 * (region.width()/region.height()) * ( (m_longlat_bbox.width()/dims.x()) / (m_longlat_bbox.height()/dims.y()) );

    double bottom_lat = m_longlat_bbox.max().y() - region.max().y()*m_longlat_bbox.height() / dims.y();
    double top_lat    = m_longlat_bbox.max().y() - region.min().y()*m_longlat_bbox.height() / dims.y();
    
    bool top_merge    = ( bottom_lat > 0 ) && ( ( 1.0 / cos(M_PI/180 * bottom_lat) ) > aspect_ratio );
    bool bottom_merge = ( top_lat    < 0 ) && ( ( 1.0 / cos(M_PI/180 * top_lat   ) ) > aspect_ratio );

    if( top_merge ) {
      children.push_back( std::make_pair( name + "4", BBox2i( region.min(), region.max() - Vector2i(0,region.height()/2) ) ) );
    }
    else {
      children.push_back( std::make_pair( name + "0", BBox2i( (region + region.min()) / 2 ) ) );
      children.push_back( std::make_pair( name + "1", BBox2i( (region + Vector2i(region.max().x(),region.min().y())) / 2 ) ) );
    }
    if( bottom_merge ) {
      children.push_back( std::make_pair( name + "5", BBox2i( region.min() + Vector2i(0,region.height()/2), region.max() ) ) );
    }
    else {
      children.push_back( std::make_pair( name + "2", BBox2i( (region + Vector2i(region.min().x(),region.max().y())) / 2 ) ) );
      children.push_back( std::make_pair( name + "3", BBox2i( (region + region.max()) / 2 ) ) );
    }
  
    return children;
  }

  void GigapanQuadTreeConfigData::metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const {
    bool root_node = ( info.name.size() == 0 );


    if ( root_node) {
      std::ostringstream json;
      fs::path file_path( info.filepath );
      fs::path json_path = fs::path(file_path).replace_extension(".json");

      json << "{" << std::endl
           << "  \"width\": "   << qtree.get_dimensions()[0] << "," << std::endl
           << "  \"height\": "  << qtree.get_dimensions()[1] << "," << std::endl
           << "  \"nlevels\": " << qtree.get_tree_levels()   << std::endl
           << "}" << std::endl;

      fs::ofstream jsonfs(json_path);
      jsonfs << json.str();
    }
  }

  // TODO: Is this actually the right function for Gigapan?
  cartography::GeoReference GigapanQuadTreeConfig::output_georef(uint32 xresolution, uint32 yresolution) {
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
} // namespace mosaic
} // namespace vw
