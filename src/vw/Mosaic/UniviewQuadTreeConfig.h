// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file UniviewQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate Uniview overlays.
///
#ifndef __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__
#define __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  class UniviewQuadTreeConfig : public QuadTreeConfig {
  public:
    UniviewQuadTreeConfig( bool terrain = false ) : m_terrain( terrain ) {}
    virtual ~UniviewQuadTreeConfig() {}

    void set_terrain( bool terrain ) { m_terrain = terrain; }

    void configure( QuadTreeGenerator &qtree ) const;
    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );
    static boost::shared_ptr<ImageResource> terrain_tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format );

  private:
    bool m_terrain;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__
