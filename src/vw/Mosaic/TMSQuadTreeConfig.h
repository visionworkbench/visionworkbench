// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file TMSQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate TMS overlays.
///
#ifndef __VW_MOSAIC_TMSQUADTREECONFIG_H__
#define __VW_MOSAIC_TMSQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/QuadTreeConfig.h>

namespace vw {
namespace mosaic {

  // This class is overkill, but is exists by analogy to others
  // like it for consistency.
  class TMSQuadTreeConfig : public QuadTreeConfig {
  public:
    virtual ~TMSQuadTreeConfig() {}
    void configure( QuadTreeGenerator& qtree ) const;
    cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0);

    // Makes paths of the form "path/name/4/6/3.jpg"
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_TMSQUADTREECONFIG_H__
