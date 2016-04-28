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


#ifndef __VW_CARTOGRAPHY_GEOREFERENCEUTILS_H__
#define __VW_CARTOGRAPHY_GEOREFERENCEUTILS_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/Cartography/Datum.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/Exception.h>
/*
// Boost
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
*/

#include <vw/Cartography/GeoReference.h>

/// \file GeoReferenceUtils.h Tools for working with GeoReferenc objects

namespace vw {
namespace cartography {



  /// Generates a new georeference which covers a sub-region of this georeference object.
  /// - Input coordinates are pixels in the corresponding image
  GeoReference crop( GeoReference const& input,
                     double upper_left_x, double upper_left_y,
                     double width=0, double height=0 );

  /// Overload of crop() that takes a bounding box object (still in pixels)
  GeoReference crop( GeoReference const& input,
                     BBox2 const& bbox );

  /// Modify the scale in the projected coordinate to pixel coordinate transform.
  /// - A larger scale increases the number of pixels.
  GeoReference resample( GeoReference const& input, double scale_x, double scale_y );
  GeoReference resample( GeoReference const& input, double scale );

  // TODO: These functions (like all other georef functions) can fail when the lon coordinates
  //       are 360 degrees off.  This should be handled in the georef functions.

  /// Given a projected coordinate in georef1, convert it to a pixel coordinate in georef2.
  Vector2 georef_point_to_georef_pixel(Vector2 const& proj_pt1,
                                       cartography::GeoReference const& georef1,
                                       cartography::GeoReference const& georef2);

  /// Given bounding box of projected coordinates in georef1, convert it to pixel coordinates in georef2.
  BBox2 georef_point_to_georef_pixel_bbox(BBox2 point_box1,
                                          cartography::GeoReference const& georef1,
                                          cartography::GeoReference const& georef2);


  /// Given a bounding box of pixels in georef1, convert it to projected coordinates in georef2.
  BBox2 georef_pixel_to_georef_point_bbox(BBox2 pixel_box1,
                                          cartography::GeoReference const& georef1,
                                          cartography::GeoReference const& georef2);



  /// The following namespace contains functions that return GeoReferences
  /// for certain well-known output styles, such as KML (and related
  /// functions involved in doing so).
  namespace output {
    namespace kml {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        Vector2 pos      = tx.forward( pixel );
        Vector2 x_vector = tx.forward( pixel+Vector2(1,0) ) - pos;
        Vector2 y_vector = tx.forward( pixel+Vector2(0,1) ) - pos;
        double degrees_per_pixel        = (std::min)( norm_2(x_vector), norm_2(y_vector) );
        double pixels_per_circumference = 360.0 / degrees_per_pixel;
        int scale_exponent = (int) ceil( log(pixels_per_circumference)/log(2.0) );
        if (scale_exponent >= 31) scale_exponent = 30;
        return 1 << scale_exponent;
      }
    } // namespace: vw::cartography::output::kml

    namespace tms {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        // It's exactly the same as the one for KML.
        return vw::cartography::output::kml::compute_resolution(tx, pixel);
      }
    } // namespace vw::cartography::output::tms
  } // namespace vw::cartography::output


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCEUTILS_H__
