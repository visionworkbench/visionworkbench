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

#include <vw/Mosaic/KMLQuadTreeConfig.h>

#include <iomanip>
#include <sstream>
#include <boost/bind.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace mosaic {

  struct KMLQuadTreeConfigData {
    BBox2 m_longlat_bbox;
    std::string m_title;
    int m_max_lod_pixels;
    int m_draw_order_offset;
    std::string m_metadata;
    mutable std::ostringstream m_root_node_tags;

    std::string kml_latlonaltbox( BBox2 const& longlat_bbox ) const;
    std::string kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox ) const;
    std::string kml_ground_overlay( std::string const& href, BBox2 const& bbox, int draw_order, int max_lod_pixels ) const;
    BBox2 pixels_to_longlat( BBox2i const& image_bbox, Vector2i const& dimensions ) const;

    std::vector<std::pair<std::string,vw::BBox2i> > branch_func( QuadTreeGenerator const&, std::string const& name, BBox2i const& region ) const;
    void metadata_func( QuadTreeGenerator const&, QuadTreeGenerator::TileInfo const& info ) const;

  public:
    KMLQuadTreeConfigData()
      : m_longlat_bbox(-180,-90,360,180),
        m_max_lod_pixels(-1),
        m_draw_order_offset(0)
    {}
  };


  vw::mosaic::KMLQuadTreeConfig::KMLQuadTreeConfig()
    : m_data( new KMLQuadTreeConfigData() )
  {}

  void vw::mosaic::KMLQuadTreeConfig::set_longlat_bbox( BBox2 const& bbox ) {
    m_data->m_longlat_bbox = bbox;
  }

  void vw::mosaic::KMLQuadTreeConfig::set_title( std::string const& title ) {
    m_data->m_title = title;
  }

  void vw::mosaic::KMLQuadTreeConfig::set_max_lod_pixels( int32 pixels ) {
    m_data->m_max_lod_pixels = pixels;
  }

  void vw::mosaic::KMLQuadTreeConfig::set_draw_order_offset( int32 offset ) {
    m_data->m_draw_order_offset = offset;
  }

  void vw::mosaic::KMLQuadTreeConfig::set_metadata( std::string const& data ) {
    m_data->m_metadata = data;
  }

  void vw::mosaic::KMLQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_crop_images( true );
    qtree.set_image_path_func( &QuadTreeGenerator::named_tiered_image_path );
    qtree.set_metadata_func( boost::bind(&KMLQuadTreeConfigData::metadata_func,m_data,_1,_2) );
    qtree.set_branch_func( boost::bind(&KMLQuadTreeConfigData::branch_func,m_data,_1,_2,_3) );
  }

  std::string KMLQuadTreeConfigData::kml_latlonaltbox( BBox2 const& longlat_bbox ) const {
    BBox2 bbox = longlat_bbox;
    bbox.crop( m_longlat_bbox );
    std::ostringstream tag;
    tag << std::setprecision(10);
    tag << "<LatLonAltBox>"
	<< "<north>" << bbox.min().y() << "</north>"
	<< "<south>" << bbox.max().y() << "</south>"
	<< "<east>" << bbox.max().x() << "</east>"
	<< "<west>" << bbox.min().x() << "</west>"
	<< "</LatLonAltBox>";
    return tag.str();
  }

  std::string KMLQuadTreeConfigData::kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox ) const {
    std::ostringstream tag;
    tag << std::setprecision(10);
    tag << "  <NetworkLink>\n"
	<< "    <name>" + name + "</name>\n"
	<< "    <Region>" + kml_latlonaltbox(longlat_bbox) + "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
	<< "    <Link><href>" + href + "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
	<< "  </NetworkLink>\n";
    return tag.str();
  }

  std::string KMLQuadTreeConfigData::kml_ground_overlay( std::string const& href, BBox2 const& bbox, int draw_order, int max_lod_pixels ) const {
    std::ostringstream tag;
    tag << std::setprecision(10);
    tag << "  <GroundOverlay>\n"
	<< "    <Region>" << kml_latlonaltbox(bbox) << "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>" << max_lod_pixels << "</maxLodPixels></Lod></Region>\n"
	<< "    <name>" << href << "</name>\n"
	<< "    <Icon><href>" << href << "</href></Icon>\n"
	<< "    " << kml_latlonaltbox(bbox) << "\n"
	<< "    <drawOrder>" << draw_order << "</drawOrder>\n"
	<< "  </GroundOverlay>\n";
    return tag.str();
  }

  BBox2 KMLQuadTreeConfigData::pixels_to_longlat( BBox2i const& image_bbox, Vector2i const& dimensions ) const {
    double width = dimensions.x(), height = dimensions.y();
    // Fractional bounding-box
    BBox2 fbb( image_bbox.min().x()/width, image_bbox.min().y()/height,
	       image_bbox.width()/width, image_bbox.height()/height );
    BBox2 bb( m_longlat_bbox.min().x()+fbb.min().x()*m_longlat_bbox.width(),
	      m_longlat_bbox.max().y()-fbb.min().y()*m_longlat_bbox.height(),
	      fbb.width()*m_longlat_bbox.width(), -fbb.height()*m_longlat_bbox.height() );
    return bb;
  }

  void KMLQuadTreeConfigData::metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const {
    fs::path file_path( info.filepath, fs::native );
    int base_len = file_path.branch_path().native_file_string().size() + 1;
    fs::path kml_path = change_extension( file_path, ".kml" );
    fs::ofstream kml( kml_path );
    kml << std::setprecision(10);

    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml><Folder>\n";
    int num_children = 0;

    bool root_node = ( info.name.size() == 0 );

    if( root_node ) {
      if( ! m_title.empty() )
	m_root_node_tags << "  <name>" << m_title << "</name>\n";
      if( ! m_metadata.empty() )
	m_root_node_tags << "  <Metadata>" << m_metadata << "</Metadata>\n";

      BBox2 bbox = (qtree.get_crop_bbox().empty()) ? m_longlat_bbox : pixels_to_longlat( qtree.get_crop_bbox(), qtree.get_dimensions() );

      double lon = (bbox.min().x()+bbox.max().x())/2;
      double lat = (bbox.min().y()+bbox.max().y())/2;
      double range = 1e5 * (bbox.width()*cos(M_PI/180*lat)-bbox.height());

      if( range > 1.2e7 ) range = 1.2e7;
      m_root_node_tags << "  <LookAt><longitude>" << lon << "</longitude><latitude>" << lat << "</latitude><range>" << range << "</range></LookAt>\n";
      m_root_node_tags << "  <Style><ListStyle><listItemType>checkHideChildren</listItemType></ListStyle></Style>\n";
    }

    std::vector<std::pair<std::string,BBox2i> > children = qtree.branches( info.name, info.region_bbox );
    for( unsigned i=0; i<children.size(); ++i ) {
      std::string kmlfile = qtree.image_path(children[i].first) + ".kml";
      if( exists( fs::path( kmlfile, fs::native ) ) ) {
	num_children++;
	kml << kml_network_link( children[i].first.substr(children[i].first.size()-1),
				 kmlfile.substr(base_len),
				 pixels_to_longlat( children[i].second, qtree.get_dimensions() ) );
      }
    }

    int max_lod = m_max_lod_pixels;
    if( num_children == 0 ) max_lod = -1;
    int draw_order = m_draw_order_offset + info.name.size();
    kml << kml_ground_overlay( file_path.leaf() + info.filetype,
			       pixels_to_longlat( info.image_bbox, qtree.get_dimensions() ),
			       draw_order, max_lod );
    
    if( root_node ) kml << m_root_node_tags.str();
    kml << "</Folder></kml>\n";
  }

  std::vector<std::pair<std::string,vw::BBox2i> > KMLQuadTreeConfigData::branch_func( QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region ) const {
    std::vector<std::pair<std::string,vw::BBox2i> > children;
    if( region.height() > qtree.get_tile_size() ) {

      double aspect_ratio = 2 * region.width() / region.height();
      double bottom_lat = m_longlat_bbox.max().y() - region.max().y()*m_longlat_bbox.height() / qtree.get_dimensions().y();
      double top_lat = m_longlat_bbox.max().y() - region.min().y()*m_longlat_bbox.height() / qtree.get_dimensions().y();
      bool top_merge = ( bottom_lat > 0 ) && ( ( 1.0 / cos(M_PI/180 * bottom_lat) ) > aspect_ratio );
      bool bottom_merge = ( top_lat < 0 ) && ( ( 1.0 / cos(M_PI/180 * top_lat) ) > aspect_ratio );

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
    }
    return children;
  }

} // namespace mosaic
} // namespace vw
