// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file GigapanQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate quadtrees for the
/// Gigapan viewer
///
#ifndef __VW_MOSAIC_GIGAPANQUADTREECONFIG_H__
#define __VW_MOSAIC_GIGAPANQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  struct GigapanQuadTreeConfigData;

  class GigapanQuadTreeConfig : public QuadTreeConfig {
    // The implementation is stored in a shared pointer so that it can
    // be safely bound to the quadtree callbacks in colsures even if
    // this config object goes out of scope.
    boost::shared_ptr<GigapanQuadTreeConfigData> m_data;

  public:
    GigapanQuadTreeConfig();
    virtual ~GigapanQuadTreeConfig() {}

    void configure( QuadTreeGenerator& qtree ) const;
    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

    // Makes paths of the form "path/name/4/6/3.jpg"
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

    // Set the extents (in degrees) of the quadtree.
    void set_longlat_bbox( BBox2 const& bbox );
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_GIGAPANQUADTREECONFIG_H__
