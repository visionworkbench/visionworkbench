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


#ifndef __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
#define __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__

#include <vw/Plate/PlateManager.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class PolarStereoPlateManager : public PlateManager<PixelT> {
  protected:
    // Transforms an input image with input georef to polar stereo
    // graphic at most ideal matching pyramid level.
    void transform_image( cartography::GeoReference const& georef,
                          ImageViewRef<PixelT>& image,
                          TransformRef& txref, int& level ) const;

  public:
    PolarStereoPlateManager(boost::shared_ptr<PlateFile> platefile):
      PlateManager<PixelT>(platefile) {}

    // Provide a georeference that represents a pyramid level
    cartography::GeoReference georeference( int level, bool north_pole,
                                            cartography::Datum const& datum ) const;
    cartography::GeoReference georeference( int level ) const;
  };

}} // end vw::platefile

#endif//__VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
