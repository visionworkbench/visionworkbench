// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file GMapQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate Google Maps overlays.
///
#ifndef __VW_MOSAIC_GMAPQUADTREECONFIG_H__
#define __VW_MOSAIC_GMAPQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  // This class is overkill, but is exists by analogy to others
  // like it for consistency.
  class GMapQuadTreeConfig : public QuadTreeConfig {
  public:
    virtual ~GMapQuadTreeConfig() {}
    void configure( QuadTreeGenerator& qtree ) const;
    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

    // Makes paths of the form "path/name/4/6/3.jpg"
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

    // Makes DiskImageResource objects with no file extension
    static boost::shared_ptr<ImageResource> tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format );

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_GMAPQUADTREECONFIG_H__
