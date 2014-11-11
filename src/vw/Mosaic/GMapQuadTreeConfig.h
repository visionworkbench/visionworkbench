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
    // - More specifically: [path]/[qtree_name]/[level]/[x_index]/[y_index]
    // - For each level, the origin is upper left, advancing right and down (a very common indexing scheme).
    // - OR IS IT?  It seems to advance up instead of down!!!!!!!
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

    // Makes DiskImageResource objects with no file extension
    static boost::shared_ptr<DstImageResource> tile_resource( QuadTreeGenerator const& qtree, QuadTreeGenerator::TileInfo const& info, ImageFormat const& format );

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_GMAPQUADTREECONFIG_H__
