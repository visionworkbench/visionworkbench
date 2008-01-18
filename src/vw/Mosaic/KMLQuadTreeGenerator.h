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

/// \file KMLQuadTreeGenerator.h
/// 
/// A subclass of ImageQuadTreeGenerator that generates multi-resolution 
/// KML overlays.
/// 
#ifndef __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
#define __VW_MOSAIC_KMLQUADTREEGENERATOR_H__

#include <iomanip>
#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  template <class PixelT>
  class KMLQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
  {
    typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

  protected:
    BBox2 total_longlat_bbox;
    int max_lod_pixels;
    int draw_order_offset;
    std::string kml_title;
    mutable std::ostringstream root_node_tags;

    std::string kml_latlonaltbox( BBox2 const& longlat_bbox ) const {
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "<LatLonAltBox>"
          << "<north>" << longlat_bbox.min().y() << "</north>"
          << "<south>" << longlat_bbox.max().y() << "</south>"
          << "<east>" << longlat_bbox.max().x() << "</east>"
          << "<west>" << longlat_bbox.min().x() << "</west>"
          << "</LatLonAltBox>";
      return tag.str();
    }

    std::string kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox ) const {
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "  <NetworkLink>\n"
          << "    <name>" + name + "</name>\n"
          << "    <Region>" + kml_latlonaltbox(longlat_bbox) + "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
          << "    <Link><href>" + href + "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
          << "  </NetworkLink>\n";
      return tag.str();
    }

    std::string kml_ground_overlay( std::string const& href, BBox2 const& bbox, BBox2 const& rbbox, int draw_order, int max_lod_pixels ) const {
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "  <GroundOverlay>\n"
          << "    <Region>" << kml_latlonaltbox(rbbox) << "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>" << max_lod_pixels << "</maxLodPixels></Lod></Region>\n"
          << "    <name>" << href << "</name>\n"
          << "    <Icon><href>" << href << "</href></Icon>\n"
          << "    " << kml_latlonaltbox(bbox) << "\n"
          << "    <drawOrder>" << draw_order+draw_order_offset << "</drawOrder>\n"
          << "  </GroundOverlay>\n";
      return tag.str();
    }

    BBox2 pixels_to_longlat( BBox2i const& image_bbox ) const {
      double width = base_type::m_source.cols(), height = base_type::m_source.rows();
      // Fractional bounding-box
      BBox2 fbb( image_bbox.min().x()/width, image_bbox.min().y()/height,
                 image_bbox.width()/width, image_bbox.height()/height );
      BBox2 bb( total_longlat_bbox.min().x()+fbb.min().x()*total_longlat_bbox.width(),
                total_longlat_bbox.max().y()-fbb.min().y()*total_longlat_bbox.height(),
                fbb.width()*total_longlat_bbox.width(), -fbb.height()*total_longlat_bbox.height() );
      return bb;
    }

  public:
    // Constructor, templatized on the source image type.  The supplied 
    // bounding box sets the size in degrees of the KML overlay, and has 
    // nothing in particular to do with the source image or pixels.
    template <class ImageT>
    KMLQuadTreeGenerator( std::string const& tree_name,
                          ImageViewBase<ImageT> const& source,
                          BBox2 total_longlat_bbox )
      : vw::mosaic::ImageQuadTreeGenerator<PixelT>( tree_name, source ),
        total_longlat_bbox( total_longlat_bbox ),
        max_lod_pixels( -1 ),
        draw_order_offset( 0 )
    {}
    
    void set_max_lod_pixels( int pixels ) {
      max_lod_pixels = pixels;
    }

    void set_draw_order_offset( int offset ) {
      draw_order_offset = offset;
    }

    void set_kml_title( std::string const& name ) {
      kml_title = name;
    }

    virtual ~KMLQuadTreeGenerator() {}
    
    virtual std::string compute_image_path( std::string const& name ) const {
      fs::path path( base_type::m_base_dir, fs::native );

      if( name.length() == 1 ) {
        path /= change_extension( path, "" ).leaf();
      }
      else {
        for ( int i=1; i<(int)name.length() - (int)base_type::m_levels_per_directory; i+=base_type::m_levels_per_directory ) {
          path /= name.substr( i, base_type::m_levels_per_directory );
        }
        path /= name.substr( 1, name.length()-1 );
      }

      return path.native_file_string();
    }

    virtual void write_meta_file( typename base_type::PatchInfo const& info ) const
    {
      bool root_node = ( info.level == base_type::m_tree_levels-1 );

      if( root_node ) {
        if( ! kml_title.empty() )
          root_node_tags << "  <name>" << kml_title << "</name>\n";
        BBox2 bbox = pixels_to_longlat( base_type::m_crop_bbox );
        double lon = (bbox.min().x()+bbox.max().x())/2;
        double lat = (bbox.min().y()+bbox.max().y())/2;
        double range = 1e5 * (bbox.width()*cos(M_PI/180*lat)-bbox.height());
        if( range > 1.2e7 ) range = 1.2e7;
        root_node_tags << "  <LookAt><longitude>" << lon << "</longitude><latitude>" << lat << "</latitude><range>" << range << "</range></LookAt>\n";
        root_node_tags << "  <Style><ListStyle><listItemType>checkHideChildren</listItemType></ListStyle></Style>\n";
      }

      BBox2 bb = pixels_to_longlat( info.image_bbox );
      BBox2 rbb = pixels_to_longlat( info.region_bbox );

      fs::path file_path( info.filepath, fs::native );
      int base_len = file_path.branch_path().native_file_string().size() + 1;

      fs::path kml_path = change_extension(file_path,".kml");
      fs::ofstream kml( kml_path );
      kml << std::setprecision(10);

      kml << "<Folder>\n";
      int children = 0;
      std::string kml0 = compute_image_path( info.name+"0" ) + ".kml";
      if( exists( fs::path( kml0, fs::native ) ) ) {
        ++children;
        kml << kml_network_link( "0", kml0.substr(base_len), (rbb+Vector2(rbb.min().x(),rbb.min().y()))/2.0 );
      }
      std::string kml1 = compute_image_path( info.name+"1" ) + ".kml";
      if( exists( fs::path( kml1, fs::native ) ) ) {
        ++children;
        kml << kml_network_link( "1", kml1.substr(base_len), (rbb+Vector2(rbb.max().x(),rbb.min().y()))/2.0 );
      }
      std::string kml2 = compute_image_path( info.name+"2" ) + ".kml";
      if( exists( fs::path( kml2, fs::native ) ) ) {
        ++children;
        kml << kml_network_link( "2", kml2.substr(base_len), (rbb+Vector2(rbb.min().x(),rbb.max().y()))/2.0 );
      }
      std::string kml3 = compute_image_path( info.name+"3" ) + ".kml";
      if( exists( fs::path( kml3, fs::native ) ) ) {
        ++children;
        kml << kml_network_link( "3", kml3.substr(base_len), (rbb+Vector2(rbb.max().x(),rbb.max().y()))/2.0 );
      }
      int max_lod = max_lod_pixels;
      if( ! children ) max_lod = -1;
      if( base_type::m_output_image_type == "auto" && extension(file_path)==".jpg" ) max_lod = -1;
      kml << kml_ground_overlay( file_path.leaf(), bb, rbb, base_type::m_tree_levels-info.level, max_lod );
      if( root_node ) kml << root_node_tags.str();
      kml << "</Folder>\n";
    }

  };


  /// A transform functor that relates unprojected lon/lat 
  /// coordinates in degrees to an unprojected pixel space 
  /// cooresponding to a standard global KML image quad-tree.
  class GlobalKMLTransform : public TransformHelper<GlobalKMLTransform,ConvexFunction,ConvexFunction>
  {
    int xresolution, yresolution;
  public:
    GlobalKMLTransform( int resolution ) : xresolution(resolution), yresolution(resolution) {}
    GlobalKMLTransform( int xresolution, int yresolution ) : xresolution(xresolution), yresolution(yresolution) {}
    
    // Convert degrees lat/lon to pixel location
    inline Vector2 forward( Vector2 const& p ) const {
      return Vector2( xresolution*(p.x()+180.0)/360.0-0.5, yresolution*(180.0-p.y())/360.0-0.5 );
    }
    
    // Convert pixel location to degrees lat/lon
    inline Vector2 reverse( Vector2 const& p ) const {
      return Vector2( 360.0*(p.x()+0.5)/xresolution-180.0, 180.0-360.0*(p.y()+0.5)/yresolution );
    }

    template <class TransformT>
    static inline int compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
      Vector2 pos = tx.forward( pixel );
      Vector2 x_vector = tx.forward( pixel+Vector2(1,0) ) - pos;
      Vector2 y_vector = tx.forward( pixel+Vector2(0,1) ) - pos;
      double degrees_per_pixel = (std::min)( norm_2(x_vector), norm_2(y_vector) );
      double pixels_per_circumference = 360.0 / degrees_per_pixel;
      int scale_exponent = (int) ceil( log(pixels_per_circumference)/log(2) );
      int resolution = 1 << scale_exponent;
      return resolution;
    }
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
