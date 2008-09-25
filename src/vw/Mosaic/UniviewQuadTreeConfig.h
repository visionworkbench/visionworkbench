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

/// \file UniviewQuadTreeConfig.h
/// 
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate Uniview overlays.
/// 
#ifndef __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__
#define __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  class UniviewQuadTreeConfig {
  public:
    UniviewQuadTreeConfig( bool terrain = false ) : m_terrain( terrain ) {}

    void set_terrain( bool terrain ) { m_terrain = terrain; }

    void configure( QuadTreeGenerator &qtree ) const;

    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );
    static boost::shared_ptr<ImageResource> terrain_tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format );

  private:
    bool m_terrain;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__
