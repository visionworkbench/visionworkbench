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


#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Cartography/GeoTransform.h>

//#include <vw/Math/BresenhamLine.h>
//#include <vw/Cartography/GeoReferenceResourcePDS.h>
//#include <vw/FileIO/DiskImageResourcePDS.h>

// Boost
//#include <boost/algorithm/string.hpp>
//#include <boost/foreach.hpp>

namespace vw {
namespace cartography {


GeoReference crop( GeoReference const& input,
                   double upper_left_x, double upper_left_y,
                   double /*width*/, double /*height*/ ) {
  Vector2 top_left_ll;
  if ( input.pixel_interpretation() == GeoReference::PixelAsArea ) {
    top_left_ll = input.pixel_to_point( Vector2(upper_left_x, upper_left_y ) - Vector2(0.5,0.5) );
  } else {
    top_left_ll = input.pixel_to_point( Vector2(upper_left_x, upper_left_y ) );
  }
  GeoReference output = input;      // Start with copy of current transform
  Matrix3x3 T = output.transform();
  T(0,2) = top_left_ll[0];          // Shift the translation to the crop region
  T(1,2) = top_left_ll[1];          //  (don't need to worry about width/height)
  output.set_transform(T);
  return output;
}

GeoReference crop( GeoReference const& input,
                   BBox2 const& bbox ) {
  // Redirect to the other georeference crop call
  return crop(input, bbox.min().x(), bbox.min().y(),
              bbox.width(), bbox.height());
}

GeoReference resample( GeoReference const& input, double scale_x, double scale_y ) {
  GeoReference output = input;
  Matrix3x3 T = output.transform();
  T(0,0) /= scale_x;
  T(1,1) /= scale_y;
  if ( input.pixel_interpretation() == GeoReference::PixelAsArea ) {
    Vector2 top_left_ll = input.pixel_to_point( -Vector2(0.5 / scale_x, 0.5 / scale_y) );
    T(0,2) = top_left_ll[0];
    T(1,2) = top_left_ll[1];
  }
  output.set_transform(T);
  return output;
}

GeoReference resample( GeoReference const& input, double scale ) {
  return resample(input, scale, scale );
}


/*
Vector2 georef_point_to_georef_pixel(Vector2 const& proj_pt1,
                                     cartography::GeoReference const& georef1,
                                     cartography::GeoReference const& georef2){
  Vector2 lonlat = georef1.point_to_lonlat(proj_pt1);
  return georef2.lonlat_to_pixel(lonlat);
}

BBox2 georef_point_to_georef_pixel_bbox(BBox2 point_box1,
                                        cartography::GeoReference const& georef1,
                                        cartography::GeoReference const& georef2){

  // Ensure we don't get incorrect results for empty boxes with strange corners.
  if (point_box1.empty())
    return BBox2();

  BBox2 out_box;

  double minx = point_box1.min().x(), maxx = point_box1.max().x();
  double miny = point_box1.min().y(), maxy = point_box1.max().y();
  double rangex = maxx-minx;
  double rangey = maxy-miny;

  // At the poles this won't be enough, more thought is needed.
  int num_steps = 100;
  for (int i = 0; i <= num_steps; i++) {
    double r = double(i)/num_steps;

    // left edge
    Vector2 P2 = Vector2(minx, miny + r*rangey);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}

    // right edge
    P2 = Vector2(maxx, miny + r*rangey);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}

    // bottom edge
    P2 = Vector2(minx + r*rangex, miny);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}

    // top edge
    P2 = Vector2(minx + r*rangex, maxy);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}
    
    // diag1
    P2 = Vector2(minx + r*rangex, miny + r*rangey);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}

    // diag2
    P2 = Vector2(maxx - r*rangex, miny + r*rangey);
    try { out_box.grow(georef_point_to_georef_pixel(P2, georef1, georef2)); }
    catch ( const std::exception & e ) {}
  }

  return grow_bbox_to_int(out_box);
}

// Given an image with georef1 and a portion of its pixels in
// pixel_box1, find the bounding box of pixel_box1 in projected
// point units for georef2.
BBox2 georef_pixel_to_georef_point_bbox(BBox2 pixel_box1,
                                  cartography::GeoReference const& georef1,
                                  cartography::GeoReference const& georef2){
  cartography::GeoTransform T(georef1, georef2);
  return T.forward_pixel_to_point_bbox(pixel_box1);
}
*/




}} // vw::cartography

#undef CHECK_PROJ_ERROR
