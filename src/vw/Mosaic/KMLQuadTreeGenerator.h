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

    BBox2 bbox;
    mutable unsigned draw_order;

    static std::string kml_latlonaltbox( BBox2 const b ) {
      std::ostringstream tag;
      tag << "<LatLonAltBox>"
          << "<north>" << b.min().y() << "</north>"
          << "<south>" << b.max().y() << "</south>"
          << "<east>" << b.max().x() << "</east>"
          << "<west>" << b.min().x() << "</west>"
          << "</LatLonAltBox>";
      return tag.str();
    }

    static std::string kml_network_link( std::string const& name, std::string const& href, BBox2 const& bbox ) {
      std::ostringstream tag;
      tag << "  <NetworkLink>\n"
          << "    <name>" + name + "</name>\n"
          << "    <Region>" + kml_latlonaltbox(bbox) + "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
          << "    <Link><href>" + href + "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
          << "  </NetworkLink>\n";
      return tag.str();
    }

    static std::string kml_ground_overlay( std::string const& name, std::string const& href, BBox2 const& bbox, int draw_order ) {
      std::ostringstream tag;
      tag << "  <GroundOverlay>\n"
          << "    <Region>" << kml_latlonaltbox(bbox) << "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
          << "    <name>" << name << "</name>\n"
          << "    <Icon><href>" << href << "</href></Icon>\n"
          << "    " << kml_latlonaltbox(bbox) << "\n"
          << "    <drawOrder>" << draw_order << "</drawOrder>\n"
          << "  </GroundOverlay>\n";
      return tag.str();
    }

  public:
    // Constructor, templatized on the source image type.  The supplied 
    // bounding box sets the size in degrees of the KML overlay, and has 
    // nothing in particular to do with the source image or pixels.
    template <class ImageT>
    KMLQuadTreeGenerator( std::string const& tree_name,
                          ImageViewBase<ImageT> const& source,
                          BBox2 bbox )
      : vw::mosaic::ImageQuadTreeGenerator<PixelT>( tree_name, source ),
        bbox( bbox ),
        draw_order(10000) {}
    
    virtual ~KMLQuadTreeGenerator() {}
    
    virtual void write_meta_file( std::string const& name,
                                  unsigned scale, 
                                  BBox2i const& image_bbox,
                                  BBox2i const& visible_bbox ) const
    {
      boost::filesystem::path file_path( name );
      std::string leaf = file_path.leaf();
      std::ofstream kml( (name+".kml").c_str() );
      float width = base_type::m_source.cols(), height = base_type::m_source.rows();
      // Fractional bounding-box
      BBox2 fbb( image_bbox.min().x()/width, image_bbox.min().y()/height,
                 image_bbox.width()/width, image_bbox.height()/height );
      BBox2 bb( bbox.min().x()+fbb.min().x()*bbox.width(),
                bbox.max().y()-fbb.min().y()*bbox.height(),
                fbb.width()*bbox.width(), -fbb.height()*bbox.height() );
      int children = 0;
      kml << "<Folder>\n";
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
      if( children > 1 ) {
        kml << kml_ground_overlay( name, leaf+"."+KMLQuadTreeGenerator::m_output_image_type, bb, draw_order-- );
      }
      else {
        remove( boost::filesystem::path( name + "." + base_type::m_output_image_type ) );
      }
      kml << "</Folder>\n";
    }

  };


} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
