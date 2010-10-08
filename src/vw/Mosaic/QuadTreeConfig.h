// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_MOSAIC_QUADTREECONFIG_H__
#define __VW_MOSAIC_QUADTREECONFIG_H__

#include <vw/Cartography/GeoReference.h>
#include <boost/shared_ptr.hpp>

namespace vw {
namespace mosaic {
  class QuadTreeGenerator;

  class QuadTreeConfig {
  public:
    virtual ~QuadTreeConfig() {}
    virtual void configure( QuadTreeGenerator& qtree ) const = 0;
    virtual cartography::GeoReference output_georef(uint32 xresolution, uint32 yresolution = 0) = 0;
    static boost::shared_ptr<QuadTreeConfig> make(const std::string& type);
  };

}} // namespace vw::mosaic

#endif
