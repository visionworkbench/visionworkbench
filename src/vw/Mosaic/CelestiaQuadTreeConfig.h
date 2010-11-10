// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CelestiaQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate Celestia "virtual textures".
///
/// Helpful doc on the subject: "Creating Textures for Celestia"
/// <http://www.lns.cornell.edu/~seb/celestia/textures.html>
///
#ifndef __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__
#define __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  // This class is overkill, but it exists by analogy to others
  // like it for consistency.
  class CelestiaQuadTreeConfig : public QuadTreeConfig {
    std::string m_module_name;
  public:
    virtual ~CelestiaQuadTreeConfig() {}
    void configure( QuadTreeGenerator& qtree ) const;
    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

    // Makes paths of the form "path/name/level1/tx_2_1.jpg"
    // 2_1 is the tile at (x,y) location (2,1), (0,0) is upper-left
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

    void metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const;

    void set_module(const std::string& module);
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__
