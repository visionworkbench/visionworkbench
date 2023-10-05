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

#ifndef __VW_CARTOGRAPHY_GEOREFERENCE_BASE_UTILS_H__
#define __VW_CARTOGRAPHY_GEOREFERENCE_BASE_UTILS_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Core/Exception.h>

/// \file GeoReferenceBaseUtils.h
/// Very small standalone utilities for working with georeferences and images

namespace vw {
namespace cartography {

vw::Vector3 projToEcef(GeoReference const& georef, vw::Vector3 const& proj);
vw::Vector3 ecefToProj(GeoReference const& georef, vw::Vector3 const& ecef);

// A function that will read a geo-referenced image, its nodata value,
// and the georeference, and will return a PixelMasked image, the nodata
// value, and the georeference.
void readGeorefImage(std::string const& image_file, 
  float & nodata_val, vw::cartography::GeoReference & georef,
  vw::ImageViewRef<vw::PixelMask<float>> & masked_image);

// Given a georeferenced image and a point in ECEF, return the closest
// pixel value in the image. If the point is outside the image, return
// an invalid pixel.
vw::PixelMask<float> closestPixelVal(vw::ImageViewRef<vw::PixelMask<float>> const& image,
                                     vw::cartography::GeoReference const& georef,
                                     vw::Vector3 const& ecef);

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_BASE_UTILS_H__
