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

#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  template <class PixelT>
  class KMLQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
  {
    typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

    BBox2 total_longlat_bbox;
    int max_lod_pixels;
    mutable unsigned draw_order;

    std::string kml_latlonaltbox( BBox2 const& longlat_bbox ) const {
      std::ostringstream tag;
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
      tag << "  <NetworkLink>\n"
          << "    <name>" + name + "</name>\n"
          << "    <Region>" + kml_latlonaltbox(longlat_bbox) + "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
          << "    <Link><href>" + href + "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
          << "  </NetworkLink>\n";
      return tag.str();
    }

    std::string kml_ground_overlay( std::string const& name, std::string const& href, BBox2 const& bbox, int draw_order, int max_lod_pixels ) const {
      std::ostringstream tag;
      tag << "  <GroundOverlay>\n"
          << "    <Region>" << kml_latlonaltbox(bbox) << "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>" << max_lod_pixels << "</maxLodPixels></Lod></Region>\n"
          << "    <name>" << name << "</name>\n"
          << "    <Icon><href>" << href << "</href></Icon>\n"
          << "    " << kml_latlonaltbox(bbox) << "\n"
          << "    <drawOrder>" << draw_order << "</drawOrder>\n"
          << "  </GroundOverlay>\n";
      return tag.str();
    }

    BBox2 pixels_to_longlat( BBox2i const& image_bbox ) const {
      float width = base_type::m_source.cols(), height = base_type::m_source.rows();
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
        draw_order( 10000 ) {}
    
    void set_max_lod_pixels( int pixels ) {
      max_lod_pixels = pixels;
    }

    virtual ~KMLQuadTreeGenerator() {}
    
    virtual void write_meta_file( std::string const& name,
                                  unsigned scale, 
                                  BBox2i const& image_bbox,
                                  BBox2i const& visible_bbox ) const
    {
      boost::filesystem::path file_path( name );
      std::string leaf = file_path.leaf();
      std::ofstream kml( (name+".kml").c_str() );
      BBox2 bb = pixels_to_longlat( image_bbox );
      kml << "<Folder>\n";
      int children = 0;
      if( exists( file_path/"0" ) ) {
        ++children;
        kml << kml_network_link( name+"/0", leaf+"/0.kml", (bb+Vector2(bb.min().x(),bb.min().y()))/2.0 );
      }
      if( exists( file_path/"1" ) ) {
        ++children;
        kml << kml_network_link( name+"/1", leaf+"/1.kml", (bb+Vector2(bb.max().x(),bb.min().y()))/2.0 );
      }
      if( exists( file_path/"2" ) ) {
        ++children;
        kml << kml_network_link( name+"/2", leaf+"/2.kml", (bb+Vector2(bb.min().x(),bb.max().y()))/2.0 );
      }
      if( exists( file_path/"3" ) ) {
        ++children;
        kml << kml_network_link( name+"/3", leaf+"/3.kml", (bb+Vector2(bb.max().x(),bb.max().y()))/2.0 );
      }
      kml << kml_ground_overlay( name, leaf+"."+KMLQuadTreeGenerator::m_output_image_type, bb, draw_order--, (children==0)?(-1):max_lod_pixels );
      if( scale == unsigned(1<<(base_type::m_tree_levels-1)) ) {
        BBox2 bbox = pixels_to_longlat( base_type::m_crop_bbox );
        double lon = (bbox.min().x()+bbox.max().x())/2;
        double lat = (bbox.min().y()+bbox.max().y())/2;
        double range = 6e4 * (bbox.width()*cos(M_PI/180*lat)-bbox.height());
        kml << "  <LookAt><longitude>" << lon << "</longitude><latitude>" << lat << "</latitude><range>" << range << "</range></LookAt>\n";
      }
      kml << "</Folder>\n";
    }

  };


} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
