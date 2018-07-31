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


#include <vw/Math/BresenhamLine.h>
#include <vw/Cartography/GeoTransform.h>

// Vision Workbench
#include <vw/Image/ImageView.h>

// Proj.4
#include <proj_api.h>

// Macro for checking Proj.4 output, something we do a lot of.
#define CHECK_PROJ_ERROR(ctx_input) if(ctx_input.error_no()) vw_throw(ProjectionErr() << "Bad projection in GeoTransform.cc. Proj.4 error: " << pj_strerrno(ctx_input.error_no()))

namespace vw {
namespace cartography {

  using vw::math::BresenhamLine;

  // Constructor
  GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                             BBox2 const& src_bbox, BBox2 const& dst_bbox) :
    m_src_georef(src_georef), m_dst_georef(dst_georef),
    m_src_bbox(src_bbox), m_dst_bbox(dst_bbox) {
    
    const std::string src_datum = m_src_georef.datum().proj4_str();
    const std::string dst_datum = m_dst_georef.datum().proj4_str();

    // This optimizes in the common case where the two images are
    // already in the same map projection, and we need only apply
    // the affine transform.
    if ((m_src_georef.overall_proj4_str() == m_dst_georef.overall_proj4_str()) &&
        (src_georef.is_lon_center_around_zero() == dst_georef.is_lon_center_around_zero()) )
      m_skip_map_projection = true;
    else
      m_skip_map_projection = false;

    // This optimizes the case where the two datums are the same,
    // and thus we don't need to call proj to convert between them
    // as we transform.
    if(src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {     
      // Set up the various variables for proj.
      m_skip_datum_conversion = false;

      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the
      // source or destination georef uses.
      ss_src << "+proj=longlat " << src_datum;
      m_src_datum_proj = ProjContext( ss_src.str() );

      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=longlat " << dst_datum;
      m_dst_datum_proj = ProjContext( ss_dst.str() );
    }
    // Because GeoTransform is typically very slow, we default to a tolerance
    // of 0.1 pixels to allow ourselves to be approximated.
    set_tolerance( 0.1 );
  }

  GeoTransform::GeoTransform(GeoTransform const& other) {
    *this = other;
  }

  GeoTransform& GeoTransform::operator=(GeoTransform const& other) {
    m_src_georef            = other.m_src_georef;
    m_dst_georef            = other.m_dst_georef;
    m_src_datum_proj        = other.m_src_datum_proj;
    m_dst_datum_proj        = other.m_dst_datum_proj;
    m_skip_map_projection   = other.m_skip_map_projection;
    m_skip_datum_conversion = other.m_skip_datum_conversion;
    return *this;
  }


  Vector2 GeoTransform::reverse(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_src_georef.point_to_pixel(m_dst_georef.pixel_to_point(v));
    Vector2 dst_lonlat = m_dst_georef.pixel_to_lonlat(v);
    if(m_skip_datum_conversion)
      return m_src_georef.lonlat_to_pixel(dst_lonlat);
    Vector2 src_lonlat = lonlat_to_lonlat(dst_lonlat, false);
    return m_src_georef.lonlat_to_pixel(src_lonlat);
  }


  Vector2 GeoTransform::pixel_to_pixel(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(m_src_georef.pixel_to_point(v));
    Vector2 src_lonlat = m_src_georef.pixel_to_lonlat(v);
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_pixel(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_pixel(dst_lonlat);
  }


  BBox2i GeoTransform::forward_bbox( BBox2i const& bbox ) const {

    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    gen_bd_and_diag_pts(bbox, points);
    
    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow( this->forward( points[ptiter] ) );
      }catch ( const std::exception & e ) {}
    }

    return grow_bbox_to_int(r);
  }

  BBox2i GeoTransform::reverse_bbox( BBox2i const& bbox ) const {
    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    gen_bd_and_diag_pts(bbox, points);

    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow( this->reverse( points[ptiter] ) );
      }catch ( const std::exception & e ) {}
    }

    return grow_bbox_to_int(r);
  }


  Vector2 GeoTransform::point_to_point( Vector2 const& v ) const {
    if (m_skip_map_projection)
      return v;
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(v);
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_point(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_point(dst_lonlat);
  }


  Vector2 GeoTransform::pixel_to_point( Vector2 const& v ) const {
    Vector2 src_point = m_src_georef.pixel_to_point(v);
    if (m_skip_map_projection)
      return src_point;
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(src_point);
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_point(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_point(dst_lonlat);
  }

  Vector2 GeoTransform::point_to_pixel( Vector2 const& v ) const {
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(v);
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(v);
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_pixel(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_pixel(dst_lonlat);
  }

  // Performs a forward or reverse datum conversion.
  Vector2 GeoTransform::lonlat_to_lonlat(Vector2 const& lonlat, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlat;

    double lon = lonlat[0] * DEG_TO_RAD; // proj4 requires radians
    double lat = lonlat[1] * DEG_TO_RAD;
    double alt = 0;

    Mutex::WriteLock write_lock(m_mutex);

    if(forward) // src to dst
      pj_transform(m_src_datum_proj.proj_ptr(), m_dst_datum_proj.proj_ptr(), 1, 0, &lon, &lat, &alt);
    else // dst to src
      pj_transform(m_dst_datum_proj.proj_ptr(), m_src_datum_proj.proj_ptr(), 1, 0, &lon, &lat, &alt);
    CHECK_PROJ_ERROR( m_src_datum_proj );
    CHECK_PROJ_ERROR( m_dst_datum_proj );

    return Vector2(lon*RAD_TO_DEG, lat*RAD_TO_DEG);
  }

  Vector3 GeoTransform::lonlatalt_to_lonlatalt(Vector3 const& lonlatalt, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlatalt;

    double lon = lonlatalt[0] * DEG_TO_RAD; // proj4 requires radians
    double lat = lonlatalt[1] * DEG_TO_RAD;
    double alt = lonlatalt[2];

    Mutex::WriteLock write_lock(m_mutex);

    if(forward) // src to dst
      pj_transform(m_src_datum_proj.proj_ptr(), m_dst_datum_proj.proj_ptr(), 1, 0, &lon, &lat, &alt);
    else // dst to src
      pj_transform(m_dst_datum_proj.proj_ptr(), m_src_datum_proj.proj_ptr(), 1, 0, &lon, &lat, &alt);
    CHECK_PROJ_ERROR( m_src_datum_proj );
    CHECK_PROJ_ERROR( m_dst_datum_proj );

    return Vector3(lon*RAD_TO_DEG, lat*RAD_TO_DEG, alt);
  }


  bool GeoTransform::check_bbox_wraparound() const {

    // Check if we are converting between georefs with different lon centers.
    bool lon_center_mismatch = (m_src_georef.is_lon_center_around_zero() != 
                                m_dst_georef.is_lon_center_around_zero()   );
    if (!lon_center_mismatch)
      return false; // No chance of wraparound with the same center.

    if (m_src_bbox.empty() || m_dst_bbox.empty())
      vw_throw(LogicErr() << "Cannot check bbox GeoTransform wraparound without both bboxes!");

    // The danger is that a small bounding box in one center can convert into
    //  a planet-spanning bounding box when converted to the other center.

    // Get the corners of the src image
    std::vector<Vector2> corners(4);
    corners[0] = m_src_bbox.min();
    corners[1] = Vector2(m_src_bbox.width()-1,0);
    corners[2] = Vector2(0, m_src_bbox.height()-1);
    corners[3] = Vector2(m_src_bbox.width()-1,m_src_bbox.height()-1);

    double center;
    if (m_src_georef.is_lon_center_around_zero()) {
      // [-180,180] to [0,360]
      center = 0;
    } else {
      // [0,360] to [-180,180]
      center = 180; 
    }
    
    bool hasLeft  = false;
    bool hasRight = false;
    for (size_t i=0; i<4; ++i) {
      Vector2 lonlat = m_src_georef.pixel_to_lonlat(corners[i]);
      double lon = lonlat[0];
      if (lon < center)
        hasLeft = true;
      else
        hasRight = true;
    }
    return (hasRight && hasLeft);
  }



  BBox2 GeoTransform::lonlat_to_lonlat_bbox( BBox2 const& bbox ) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (bbox.empty())
      return BBox2();

    BBox2 out_box;
    
    double minx = bbox.min().x(), maxx = bbox.max().x();
    double miny = bbox.min().y(), maxy = bbox.max().y();
    double rangex = maxx-minx;
    double rangey = maxy-miny;

    // At the poles this won't be enough, more thought is needed.
    int num_steps = 100;
    for (int i = 0; i <= num_steps; i++) {
      double r = double(i)/num_steps;

      // left edge
      Vector2 P2 = Vector2(minx, miny + r*rangey);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}

      // right edge
      P2 = Vector2(maxx, miny + r*rangey);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}

      // bottom edge
      P2 = Vector2(minx + r*rangex, miny);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}

      // top edge
      P2 = Vector2(minx + r*rangex, maxy);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}
      
      // diag1
      P2 = Vector2(minx + r*rangex, miny + r*rangey);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}

      // diag2
      P2 = Vector2(maxx - r*rangex, miny + r*rangey);
      try { out_box.grow(lonlat_to_lonlat(P2)); }
      catch ( const std::exception & e ) {}
    }

    return out_box;
  }
  
  BBox2 GeoTransform::pixel_to_point_bbox( BBox2 const& pixel_bbox ) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (pixel_bbox.empty())
      return BBox2();
    
    std::vector<vw::Vector2> points;
    gen_bd_and_diag_pts(pixel_bbox, points);
 
    BBox2 point_bbox;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        point_bbox.grow(this->pixel_to_point(points[ptiter]));
      }catch ( const std::exception & e ) {}
    }
    
    return point_bbox;
  }

  BBox2 GeoTransform::point_to_pixel_bbox( BBox2 const& point_bbox ) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (point_bbox.empty())
      return BBox2();

    BBox2 out_box;

    double minx = point_bbox.min().x(), maxx = point_bbox.max().x();
    double miny = point_bbox.min().y(), maxy = point_bbox.max().y();
    double rangex = maxx-minx;
    double rangey = maxy-miny;

    // At the poles this won't be enough, more thought is needed.
    int num_steps = 100;
    for (int i = 0; i <= num_steps; i++) {
      double r = double(i)/num_steps;

      // left edge
      Vector2 P2 = Vector2(minx, miny + r*rangey);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}

      // right edge
      P2 = Vector2(maxx, miny + r*rangey);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}

      // bottom edge
      P2 = Vector2(minx + r*rangex, miny);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}

      // top edge
      P2 = Vector2(minx + r*rangex, maxy);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}
      
      // diag1
      P2 = Vector2(minx + r*rangex, miny + r*rangey);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}

      // diag2
      P2 = Vector2(maxx - r*rangex, miny + r*rangey);
      try { out_box.grow(point_to_pixel(P2)); }
      catch ( const std::exception & e ) {}
    }

    return grow_bbox_to_int(out_box);
  }

  std::ostream& operator<<(std::ostream& os, const GeoTransform& trans) {   
    os << "Source Georef  : \n"       << trans.m_src_georef << "\n";
    os << "Dest Georef  : \n"         << trans.m_dst_georef << "\n";
    os << "skip_map_projection: "     << trans.m_skip_map_projection << "\n";
    os << "skip_datum_conversion: " << trans.m_skip_datum_conversion << "\n";
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
}} // namespace vw::cartography

#undef CHECK_PROJ_ERROR
