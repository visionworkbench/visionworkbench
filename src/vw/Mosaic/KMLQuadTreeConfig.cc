// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Mosaic/KMLQuadTreeConfig.h>

#include <vw/Image.h>
#include <vw/FileIO/DiskImageResourcePNG.h>

#include <iomanip>
#include <sstream>
#include <boost/bind.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using vw::cartography::GeoReference;

namespace vw {

  // This wrapper class intercepts premultiplied alpha data being written
  // to a PNG resource, and implements a pyramid-based hole-filling
  // algorithm to extrapolate data into the alpha-masked regions of the
  // image.
  //
  // This is a workaround hack for a Google Earth bug, in which GE's
  // rendering of semi-transparent GroundOverlays interpolates
  // alpha-masked (i.e. invalid) data, resulting in annoying (generally
  // black) fringes around semi-transparent images.
  class DiskImageResourcePNGAlphaHack : public DiskImageResourcePNG {
  public:

    DiskImageResourcePNGAlphaHack( std::string const& filename, ImageFormat const& format ) : DiskImageResourcePNG(filename,format) {}

    void write( ImageBuffer const& src, BBox2i const& bbox ) {
      int levels = (int) floor(((std::min)(log((double)bbox.width()),log((double)bbox.height())))/log(2.));
      if( levels<2 || src.unpremultiplied || !(src.format.pixel_format==VW_PIXEL_RGBA || src.format.pixel_format==VW_PIXEL_GRAYA) )
        return DiskImageResourcePNG::write(src,bbox);

      std::vector<ImageView<PixelRGBA<float> > > pyramid(levels);
      pyramid[0].set_size( bbox.width(), bbox.height() );
      convert( pyramid[0].buffer(), src, m_rescale );

      std::vector<float> kernel(2);
      kernel[0] = kernel[1] = 0.5;
      for( int i=1; i<levels; ++i ) {
        pyramid[i] = subsample(separable_convolution_filter(pyramid[i-1],kernel,kernel,1,1,ConstantEdgeExtension()),2);
      }

      for( int i=levels-1; i>0; --i ) {
        ImageView<PixelRGBA<float> > up = resample(pyramid[i],2.0,pyramid[i-1].cols(),pyramid[i-1].rows(),ConstantEdgeExtension(),BilinearInterpolation());
        if(i>1) {
          for( int y=0; y<up.rows(); ++y ) {
            for( int x=0; x<up.cols(); ++x ) {
              if(pyramid[i-1](x,y).a()==0.0) {
                pyramid[i-1](x,y) = up(x,y);
              }
            }
          }
        }
        else {
          for( int y=0; y<up.rows(); ++y ) {
            for( int x=0; x<up.cols(); ++x ) {
              if(pyramid[0](x,y).a()==0.0) {
                if(up(x,y).a()!=0.0) {
                  pyramid[0](x,y) = up(x,y) / up(x,y).a();
                  pyramid[0](x,y).a() = 0;
                }
              }
              else {
                pyramid[0](x,y).r() /= pyramid[0](x,y).a();
                pyramid[0](x,y).g() /= pyramid[0](x,y).a();
              pyramid[0](x,y).b() /= pyramid[0](x,y).a();
              }
            }
          }
        }
      }

      ImageBuffer buffer = pyramid[0].buffer();
      buffer.unpremultiplied = true;
      DiskImageResourcePNG::write(buffer,bbox);
    }
  };

namespace mosaic {

  struct KMLQuadTreeConfigData {
    BBox2 m_longlat_bbox;
    std::string m_title;
    int m_max_lod_pixels;
    int m_draw_order_offset;
    std::string m_metadata;
    mutable std::ostringstream m_root_node_tags;

    std::string kml_latlonbox( BBox2 const& longlat_bbox, bool alt ) const;
    std::string kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox, int min_lod_pixels ) const;
    std::string kml_ground_overlay( std::string const& href, BBox2 const& region_bbox, BBox2 const& image_bbox, int draw_order, int min_lod_pixels, int max_lod_pixels ) const;
    BBox2 pixels_to_longlat( BBox2i const& image_bbox, Vector2i const& dimensions ) const;

    std::vector<std::pair<std::string,vw::BBox2i> > branch_func( QuadTreeGenerator const&, std::string const& name, BBox2i const& region ) const;
    void metadata_func( QuadTreeGenerator const&, QuadTreeGenerator::TileInfo const& info ) const;
    boost::shared_ptr<DstImageResource> tile_resource_func( QuadTreeGenerator const&, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) const;

  public:
    KMLQuadTreeConfigData()
      : m_longlat_bbox(-180,-90,360,180),
        m_max_lod_pixels(1024),
        m_draw_order_offset(0)
    {}
  };


  KMLQuadTreeConfig::KMLQuadTreeConfig()
    : m_data( new KMLQuadTreeConfigData() )
  {}

  void KMLQuadTreeConfig::set_longlat_bbox( BBox2 const& bbox ) {
    m_data->m_longlat_bbox = bbox;
  }

  void KMLQuadTreeConfig::set_title( std::string const& title ) {
    m_data->m_title = title;
  }

  void KMLQuadTreeConfig::set_max_lod_pixels( int32 pixels ) {
    m_data->m_max_lod_pixels = pixels;
  }

  void KMLQuadTreeConfig::set_draw_order_offset( int32 offset ) {
    m_data->m_draw_order_offset = offset;
  }

  void KMLQuadTreeConfig::set_metadata( std::string const& data ) {
    m_data->m_metadata = data;
  }

  void KMLQuadTreeConfig::configure( QuadTreeGenerator& qtree ) const {
    qtree.set_cull_images( true );
    qtree.set_file_type( "auto" );
    qtree.set_image_path_func( QuadTreeGenerator::named_tiered_image_path() );
    qtree.set_metadata_func( boost::bind(&KMLQuadTreeConfigData::metadata_func,m_data,_1,_2) );
    qtree.set_branch_func( boost::bind(&KMLQuadTreeConfigData::branch_func,m_data,_1,_2,_3) );
    qtree.set_tile_resource_func( boost::bind(&KMLQuadTreeConfigData::tile_resource_func,m_data,_1,_2,_3) );
  }

  GeoReference KMLQuadTreeConfig::output_georef(uint32 xresolution, uint32 yresolution) {
    GeoReference r;
    r.set_pixel_interpretation(GeoReference::PixelAsArea);

    // Note: the global KML pixel space extends to +/- 180 degrees
    // latitude as well as longitude.
    Matrix3x3 transform;
    transform(0,0) = 360.0 / xresolution;
    transform(0,2) = -180;
    transform(1,1) = -360.0 / yresolution;
    transform(1,2) = 180;
    transform(2,2) = 1;
    r.set_transform(transform);

    return r;
  }

  std::string KMLQuadTreeConfigData::kml_latlonbox( BBox2 const& longlat_bbox, bool alt ) const {
    BBox2 bbox = longlat_bbox;
    std::ostringstream tag;
    tag << std::setprecision(10);

    std::string elt;
    if (alt)
      elt = "LatLonAltBox";
    else
      elt = "LatLonBox";

    tag << "<" << elt << ">"
        << "<north>" << bbox.min().y() << "</north>"
        << "<south>" << bbox.max().y() << "</south>"
        << "<east>" << bbox.max().x() << "</east>"
        << "<west>" << bbox.min().x() << "</west>"
        << "</" << elt << ">";
    return tag.str();
  }

  std::string KMLQuadTreeConfigData::kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox, int min_lod_pixels ) const {
    std::ostringstream tag;
    tag << std::setprecision(10);
    tag << "  <NetworkLink>\n"
        << "    <name>" << name << "</name>\n"
        << "    <Region>" << kml_latlonbox(longlat_bbox, true) << "<Lod><minLodPixels>" << min_lod_pixels << "</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
        << "    <Link><href>" << href << "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
        << "  </NetworkLink>\n";
    return tag.str();
  }

  std::string KMLQuadTreeConfigData::kml_ground_overlay( std::string const& href, BBox2 const& region_bbox, BBox2 const& image_bbox, int draw_order, int min_lod_pixels, int max_lod_pixels ) const {
    std::ostringstream tag;
    tag << std::setprecision(10);
    tag << "  <GroundOverlay>\n"
        << "    <Region>" << kml_latlonbox(region_bbox, true) << "<Lod><minLodPixels>" << min_lod_pixels << "</minLodPixels><maxLodPixels>" << max_lod_pixels << "</maxLodPixels></Lod></Region>\n"
        << "    <name>" << href << "</name>\n"
        << "    <Icon><href>" << href << "</href></Icon>\n"
        << "    " << kml_latlonbox(image_bbox, false) << "\n"
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
    std::ostringstream kml;
    fs::path file_path( info.filepath, fs::native );
    size_t base_len = file_path.branch_path().native_file_string().size() + 1;
    fs::path kml_path = change_extension( file_path, ".kml" );
    kml << std::setprecision(10);

    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n";
    kml << "<Folder>\n";
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
                                 pixels_to_longlat( children[i].second, qtree.get_dimensions() ),
                                 qtree.get_tile_size()/2);
      }
    }

    int max_lod = m_max_lod_pixels;
    if( num_children == 0 ) max_lod = -1;
    int draw_order = m_draw_order_offset + int(info.name.size());
    BBox2i go_bbox = (qtree.get_crop_images() ? info.image_bbox : info.region_bbox);
    if( exists( fs::path( info.filepath + info.filetype ) ) ) {
      kml << kml_ground_overlay( file_path.leaf() + info.filetype,
                                 pixels_to_longlat( info.region_bbox, qtree.get_dimensions() ),
                                 pixels_to_longlat( go_bbox, qtree.get_dimensions() ),
                                 draw_order, qtree.get_tile_size()/2, max_lod );
      num_children++;
    }

    if( root_node ) {
      kml << m_root_node_tags.str();
      num_children++;
    }

    kml << "</Folder>\n</kml>\n";

    // Skip this file if there aren't any real contents
    if( num_children != 0 ) {
      fs::ofstream kmlfs( kml_path );
      kmlfs << kml.str();
    }
  }

  std::vector<std::pair<std::string,vw::BBox2i> > KMLQuadTreeConfigData::branch_func( QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region ) const {
    std::vector<std::pair<std::string,vw::BBox2i> > children;
    if( region.height() > qtree.get_tile_size() ) {

      Vector2i dims = qtree.get_dimensions();
      double aspect_ratio = 2 * (region.width()/region.height()) * ( (m_longlat_bbox.width()/dims.x()) / (m_longlat_bbox.height()/dims.y()) );

      double bottom_lat = m_longlat_bbox.max().y() - region.max().y()*m_longlat_bbox.height() / dims.y();
      double top_lat = m_longlat_bbox.max().y() - region.min().y()*m_longlat_bbox.height() / dims.y();
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

  boost::shared_ptr<DstImageResource> KMLQuadTreeConfigData::tile_resource_func( QuadTreeGenerator const&, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format ) const {
    create_directories( fs::path( info.filepath, fs::native ).branch_path() );
    if( info.filetype == ".png" && (format.pixel_format==VW_PIXEL_RGBA || format.pixel_format==VW_PIXEL_GRAYA) ) {
      return boost::shared_ptr<DstImageResource>( new DiskImageResourcePNGAlphaHack( info.filepath+info.filetype, format ) );
    }
    else {
      return boost::shared_ptr<DstImageResource>( DiskImageResource::create( info.filepath+info.filetype, format ) );
    }
  }

} // namespace mosaic
} // namespace vw
