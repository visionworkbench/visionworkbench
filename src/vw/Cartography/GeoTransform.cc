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

// TODO(oalexan1): Wipe all mention of proj4_str.
// Everything must be done by setting the OGRSpatialReference object.

#include <vw/Math/BresenhamLine.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Image/ImageView.h>

#include <boost/algorithm/string.hpp>

namespace vw {

  // Constructor
  GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                             BBox2 const& src_bbox, BBox2 const& dst_bbox) :
    m_src_georef(src_georef), m_dst_georef(dst_georef),
    m_src_bbox(src_bbox), m_dst_bbox(dst_bbox) {

    // Coordinates pass between the georefs through lon-lat, with no datum
    // conversion, so the datums must agree.
    Datum const& src_datum = m_src_georef.datum();
    Datum const& dst_datum = m_dst_georef.datum();
    if (!same_datum(src_datum, dst_datum))
      vw_throw(ArgumentErr() << "Cannot transform between georeferences with "
               << "different datums: " << src_datum.name() << " and "
               << dst_datum.name() << ". Use gdalwarp to convert one to the "
               << "datum of the other.\n");

    // If the lonlat to lonlat transform does not change the results, skip it.
    // This results in a very notable speedup.
    if (m_src_georef.get_wkt() == m_dst_georef.get_wkt() &&
        !m_src_georef.image_ll_box().empty() &&
        !m_dst_georef.image_ll_box().empty())
      m_skip_map_projection = true;
    else
      m_skip_map_projection = false;
    
    // Watch for 360 degree offsets
    if (m_skip_map_projection) {
      vw::Vector2 src_ll = m_src_georef.image_ll_box().min();
      vw::Vector2 dst_ll = m_dst_georef.point_to_lonlat(m_dst_georef.lonlat_to_point(src_ll));
      if (norm_2(src_ll - dst_ll) > 1e-8) 
        m_skip_map_projection = false;
    }
  }

  // Try to make lon be in the range [-180,180]
  void wrapLon(double &lon) {
    if (lon > 180)
      lon -= 360;
    if (lon < -180)
      lon += 360;
  }
  
  // Inverse of forward()
  Vector2 GeoTransform::reverse(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_src_georef.point_to_pixel(m_dst_georef.pixel_to_point(v));
    Vector2 dst_lonlat = m_dst_georef.pixel_to_lonlat(v);
    
    // For a georef with projection, must ensure input lonlat is in the
    // [-180,180] range.
    if (m_src_georef.is_projected())
      wrapLon(dst_lonlat[0]);

    return m_src_georef.lonlat_to_pixel(dst_lonlat);
  }

  /// Given a pixel coordinate of an image in a source
  /// georeference frame, this routine computes the corresponding
  /// pixel in the destination (transformed) image.
  Vector2 GeoTransform::forward(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(m_src_georef.pixel_to_point(v));

    Vector2 src_lonlat = m_src_georef.pixel_to_lonlat(v);
    
    // For a geo ref with projection, must ensure input lonlat is in the
    // [-180,180] range.
    if (m_dst_georef.is_projected())
      wrapLon(src_lonlat[0]);

    return m_dst_georef.lonlat_to_pixel(src_lonlat);
  }
  
  // Apply the forward transform to a box
  BBox2i GeoTransform::forward_bbox(BBox2i const& bbox) const {

    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    sample_int_box(bbox, points);
    
    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow(this->forward(points[ptiter]));
      }catch (const std::exception & e) {}
    }

    return grow_bbox_to_int(r);
  }

  // Apply the reverse transform to a box
  BBox2i GeoTransform::reverse_bbox(BBox2i const& bbox) const {
    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    sample_int_box(bbox, points);

    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow(this->reverse(points[ptiter]));
      } catch (const std::exception & e) {}
    }
    
    return grow_bbox_to_int(r);
  }

  Vector2 GeoTransform::point_to_point(Vector2 const& v) const {
    if (m_skip_map_projection)
      return v;
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(v);
    return m_dst_georef.lonlat_to_point(src_lonlat);
  }


  Vector2 GeoTransform::pixel_to_point(Vector2 const& v) const {
    Vector2 src_point = m_src_georef.pixel_to_point(v);
    if (m_skip_map_projection)
      return src_point;
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(src_point);
    return m_dst_georef.lonlat_to_point(src_lonlat);
  }

  Vector2 GeoTransform::point_to_pixel(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(v);
     
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(v);
    return m_dst_georef.lonlat_to_pixel(src_lonlat);
  }

  BBox2 GeoTransform::point_to_point_bbox(BBox2 const& point_bbox) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (point_bbox.empty())
      return BBox2();

    std::vector<vw::Vector2> points;
    sample_float_box(point_bbox, points);

    BBox2 out_box;

    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        out_box.grow(point_to_point(points[ptiter])); 
      }catch (const std::exception & e) {}
    }
    
    return out_box;
  }

  BBox2 GeoTransform::pixel_to_point_bbox(BBox2 const& pixel_bbox) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (pixel_bbox.empty())
      return BBox2();
    
    std::vector<vw::Vector2> points;
    sample_int_box(pixel_bbox, points);
 
    BBox2 point_bbox;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        point_bbox.grow(this->pixel_to_point(points[ptiter]));
      }catch (const std::exception & e) {}
    }
    
    return point_bbox;
  }

  BBox2 GeoTransform::point_to_pixel_bbox(BBox2 const& point_bbox) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (point_bbox.empty())
      return BBox2();

    std::vector<vw::Vector2> points;
    sample_float_box(point_bbox, points);

    BBox2 out_box;

    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        out_box.grow(point_to_pixel(points[ptiter]));
      }catch (const std::exception & e) {}
    }
    
    return grow_bbox_to_int(out_box);
  }

  std::ostream& operator<<(std::ostream& os, const GeoTransform& trans) {   
    os << "Source Georef  : \n"       << trans.m_src_georef << "\n";
    os << "Dest Georef  : \n"         << trans.m_dst_georef << "\n";
    os << "skip_map_projection: "     << trans.m_skip_map_projection << "\n";
    return os;
  }

  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef) {

    GeoTransform gtx(src_georef, dst_georef);

    // Iterate over the image, transforming the first two coordinates
    // in the Vector one at a time.  The third coordinate is taken to
    // be the altitude value, and this value is not touched.
    for (int32 j=0; j < point_image.rows(); ++j) {
      for (int32 i=0; i < point_image.cols(); ++i) {
        if (point_image(i,j) != Vector3()) {
          Vector2 in(point_image(i,j)[0], point_image(i,j)[1]);
          Vector2 out = gtx.forward(in);
          point_image(i,j).x() = out[0];
          point_image(i,j).y() = out[1];
        }
      }
    }
  }
} // namespace vw
