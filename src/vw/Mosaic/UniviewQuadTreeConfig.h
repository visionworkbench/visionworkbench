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
    static boost::shared_ptr<DstImageResource> terrain_tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format );

    void metadata_func( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info ) const;
    void set_module(const std::string& module);

  private:
    bool m_terrain;
    std::string m_module_name;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_UNIVIEWQUADTREECONFIG_H__
