// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
