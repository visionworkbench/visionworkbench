// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file KMLQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate multi-resolution
/// KML overlays.
///
#ifndef __VW_MOSAIC_KMLQUADTREECONFIG_H__
#define __VW_MOSAIC_KMLQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  struct KMLQuadTreeConfigData;

  cartography::GeoReference kml_output_georef(uint32 xresolution, uint32 yresolution = 0);

  class KMLQuadTreeConfig : public QuadTreeConfig {
  public:
    KMLQuadTreeConfig();
    virtual ~KMLQuadTreeConfig() {}

    // Set the extents (in degrees) of the quadtree.
    void set_longlat_bbox( BBox2 const& bbox );

    // Set the title of the root KML file.
    void set_title( std::string const& title );

    // Set the maximum level of detail (in pixels) at which each resolution
    // of the quadtree is displayed.
    void set_max_lod_pixels( int32 pixels );

    // Set the drawOrder offset.  Overlays with positive offets are drawn on top.
    void set_draw_order_offset( int32 offset );

    // Set an option string of additional metadata to be included in the root KML file.
    void set_metadata( std::string const& data );

    // Configure the given quadtree to generate this KML.  This enables image
    // cropping and sets the image path function, branch function, and metadata
    // function.  If you intend to override any of these, be sure to do so
    // *after* calling configure() or your changes will be overwritten.
    void configure( QuadTreeGenerator& qtree ) const;

    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

  private:
    // The implementation is stored in a shared pointer so that it can
    // be safely bound to the quadtree callbacks in colsures even if
    // this config object goes out of scope.
    boost::shared_ptr<KMLQuadTreeConfigData> m_data;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREECONFIG_H__
