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

// Proj
#include <proj.h>

// Macro for checking Proj.4 output, something we do a lot of.
#define CHECK_PROJ_ERROR(ctx_input) if(ctx_input.error_no()) vw_throw(ProjectionErr() << "Bad projection in GeoTransform.cc. Proj.4 error: " << (ctx_input.error_no()))

namespace vw {
namespace cartography {

  using vw::math::BresenhamLine;

  // Constructor
  GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                             BBox2 const& src_bbox, BBox2 const& dst_bbox) :
    m_src_georef(src_georef), m_dst_georef(dst_georef),
    m_src_bbox(src_bbox), m_dst_bbox(dst_bbox) {

    // TODO(oalexan1): Must complete these given what we know about input images.
    PJ_AREA *area = NULL;
    
    // TODO(oalexan1): Need to think about threading here
    m_pj_context = proj_context_create();
    m_pj_transform = proj_create_crs_to_crs(m_pj_context,
                                            m_src_georef.overall_proj4_str().c_str(),
                                            m_dst_georef.overall_proj4_str().c_str(),
                                            area);
    if (m_pj_transform == NULL)
      vw::vw_throw(vw::ArgumentErr() << "Failed to initialize a geotransform for proj4 strings: "
                   << m_src_datum_proj.m_proj4_str << " and " 
                   << m_dst_datum_proj.m_proj4_str << ".\n");
 
//     // TODO(oalexan1): This code is repeated in 3 places
//     // https://proj-tmp.readthedocs.io/en/docs/development/quickstart.html
//     PJ* P_for_GIS = proj_normalize_for_visualization(m_pj_context, m_pj_transform);
//     if (0 == P_for_GIS) 
//       vw::vw_throw(vw::ArgumentErr() << "Failed to normalize a geotransform for proj4 strings: "
//                    << m_src_datum_proj.m_proj4_str << " and " 
//                    << m_dst_datum_proj.m_proj4_str << ".\n");
    
//     proj_destroy(m_pj_transform);
//     m_pj_transform = P_for_GIS;    
    
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
    if (src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {
      
      // Set up the various variables for proj.
      m_skip_datum_conversion = false;

      // TODO(oalexan1): Review this code, looks suspicious.
      
      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the
      // source or destination georef uses.
      ss_src << "+proj=longlat " << src_datum;
      m_src_datum_proj = ProjContext( ss_src.str() );
      
      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=longlat " << dst_datum;
      m_dst_datum_proj = ProjContext(ss_dst.str());
    }
    
    // Because GeoTransform is typically very slow, we default to a tolerance
    // of 0.1 pixels to allow ourselves to be approximated.
    set_tolerance(0.1);
    
    if (!m_skip_datum_conversion) {
      // This is a bugfix for the situation when the datums are in
      // fact the same, even if the precise datum strings differ. If
      // the change-of-datum transform returns the identity, avoid
      // doing it, as it is very slow. Here pick some samples in
      // radians, with both positive and negative longitude. Don't
      // only pick multiples of PI.

      // This check must happen after the tolerance is set
      Mutex::WriteLock write_lock(m_mutex);
      double alt = 0.0;
      double max_err = 0.0;
      for (double lon = -360; lon <= 360; lon += 90) {
        for (double lat = -90; lat <= 90; lat += 45) {
          Vector2 lon_lat(lon, lat);

          bool forward = true;
          Vector2 lon_lat2 = GeoTransform::lonlat_to_lonlat(lon_lat, forward);
          lon_lat2 = GeoTransform::lonlat_to_lonlat(lon_lat2, !forward);
          max_err = std::max(max_err, norm_2(lon_lat - lon_lat2));
        }
      }
      
      if (max_err < 1.0e-10)
        m_skip_datum_conversion = true;
    }
    
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

    // TODO(oalexan1): Must complete these given what we know about input images.
    PJ_AREA *area = NULL;
    
    // TODO(oalexan1): Need to think about threading here
    m_pj_context = proj_context_create();
    m_pj_transform = proj_create_crs_to_crs(m_pj_context,
                                            m_src_georef.overall_proj4_str().c_str(),
                                            m_dst_georef.overall_proj4_str().c_str(),
                                            area);
    if (m_pj_transform == NULL)
      vw::vw_throw(vw::ArgumentErr() << "Failed to initialize a geotransform for proj4 strings: "
                   << m_src_datum_proj.m_proj4_str << " and " 
                   << m_dst_datum_proj.m_proj4_str << ".\n");
 
//     // TODO(oalexan1): This code is repeated in 3 places
//     // https://proj-tmp.readthedocs.io/en/docs/development/quickstart.html
//     PJ* P_for_GIS = proj_normalize_for_visualization(m_pj_context, m_pj_transform);
//     if (0 == P_for_GIS) 
//       vw::vw_throw(vw::ArgumentErr() << "Failed to normalize a geotransform for proj4 strings: "
//                    << m_src_datum_proj.m_proj4_str << " and " 
//                    << m_dst_datum_proj.m_proj4_str << ".\n");
    
//     proj_destroy(m_pj_transform);
//     m_pj_transform = P_for_GIS;    
    
    return *this;
  }

  GeoTransform::~GeoTransform() {
    // TODO(oalexan1): Must figure out deallocation!
    //proj_context_destroy(m_pj_context);
    // proj_destroy(m_pj_transform);
  }

  // TODO(oalexan1): Must add a forward() function!
  
  // TODO(oalexan1): Must use directly the Proj functions,
  // without going to lon-lat!
  Vector2 GeoTransform::reverse(Vector2 const& v) const {
    Vector2 proj_pt = m_dst_georef.pixel_to_point(v);

    PJ_COORD c_in, c_out;
    c_in.lpzt.z = 0.0;
    c_in.lpzt.t = HUGE_VAL;
    c_in.enu.e = proj_pt[0]; // easting
    c_in.enu.n = proj_pt[1]; // northing

    // http://even.rouault.free.fr/proj_cpp_api/rst_generated/html/development/quickstart.html
    c_out = proj_trans(m_pj_transform, PJ_INV, c_in);

    return m_src_georef.point_to_pixel(Vector2(c_out.enu.e, c_out.enu.n));

    // TODO(oalexan1): Review here.
//     bool forward = false;
//     if (m_skip_map_projection)
//       return m_src_georef.point_to_pixel(m_dst_georef.pixel_to_point(v));
//     Vector2 dst_lonlat = m_dst_georef.pixel_to_lonlat(v);
//     if (m_skip_datum_conversion)
//       return m_src_georef.lonlat_to_pixel(dst_lonlat);
//     Vector2 src_lonlat = lonlat_to_lonlat(dst_lonlat, forward);
//     return m_src_georef.lonlat_to_pixel(src_lonlat);
  }


  // TODO(oalexan1): Must use directly the Proj functions,
  // without going to lon-lat!
  /// Given a pixel coordinate of an image in a source
  /// georeference frame, this routine computes the corresponding
  /// pixel in the destination (transformed) image.
  Vector2 GeoTransform::forward(Vector2 const& v) const {

    Vector2 proj_pt = m_src_georef.pixel_to_point(v);

    PJ_COORD c_in, c_out;
    c_in.lpzt.z = 0.0;
    c_in.lpzt.t = HUGE_VAL;
    c_in.enu.e = proj_pt[0]; // easting
    c_in.enu.n = proj_pt[1]; // northing

    // http://even.rouault.free.fr/proj_cpp_api/rst_generated/html/development/quickstart.html
    c_out = proj_trans(m_pj_transform, PJ_FWD, c_in);
    
    return m_dst_georef.point_to_pixel(Vector2(c_out.enu.e, c_out.enu.n));

    // TODO(oalexan1): Review here.
//     bool forward = true;
//     if (m_skip_map_projection)
//       return m_dst_georef.point_to_pixel(m_src_georef.pixel_to_point(v));
//     Vector2 src_lonlat = m_src_georef.pixel_to_lonlat(v);
//     if(m_skip_datum_conversion)
//       return m_dst_georef.lonlat_to_pixel(src_lonlat);
//     Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, forward);
//     return m_dst_georef.lonlat_to_pixel(dst_lonlat);
  }
  
  // TODO(oalexan1): GDAL has a function for this which is likely better written.
  BBox2i GeoTransform::forward_bbox( BBox2i const& bbox ) const {

    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    sample_int_box(bbox, points);
    
    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow( this->forward( points[ptiter] ) );
      }catch ( const std::exception & e ) {}
    }

    return grow_bbox_to_int(r);
  }

  // TODO(oalexan1): GDAL has a function for this which is likely better written.
  BBox2i GeoTransform::reverse_bbox( BBox2i const& bbox ) const {
    if (bbox.empty()) return BBox2();

    std::vector<vw::Vector2> points;
    sample_int_box(bbox, points);

    BBox2 r;
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        r.grow( this->reverse( points[ptiter] ) );
      }catch ( const std::exception & e ) {}
    }

    return grow_bbox_to_int(r);
  }


  Vector2 GeoTransform::point_to_point(Vector2 const& v) const {
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
  // TODO(oalexan1): This is a very bad function. It does not take into account
  // that when the datums differ input 0 height may result in non-zero output height.
  Vector2 GeoTransform::lonlat_to_lonlat(Vector2 const& lonlat, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlat;

    Vector3 lonlatalt(lonlat[0], lonlat[1], 0.0);
    return subvector(GeoTransform::lonlatalt_to_lonlatalt(lonlatalt, forward), 0, 2);
  }

  Vector3 GeoTransform::lonlatalt_to_lonlatalt(Vector3 const& lonlatalt, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlatalt;

    // TODO(oalexan1): Wipe this code
#if 0
    double lon = lonlatalt[0] * DEG_TO_RAD; // proj4 requires radians
    double lat = lonlatalt[1] * DEG_TO_RAD;
    double alt = lonlatalt[2];

    Mutex::WriteLock write_lock(m_mutex);
    // TODO(oalexan1): Review here.
    if(forward) // src to dst
      pj_transform(m_src_datum_proj.proj_ptr(), m_dst_datum_proj.proj_ptr(),
                   1, 0, &lon, &lat, &alt);
    else // dst to src
      pj_transform(m_dst_datum_proj.proj_ptr(), m_src_datum_proj.proj_ptr(),
                   1, 0, &lon, &lat, &alt);
    CHECK_PROJ_ERROR(m_src_datum_proj);
    CHECK_PROJ_ERROR(m_dst_datum_proj);

    return Vector3(lon*RAD_TO_DEG, lat*RAD_TO_DEG, alt);
#endif

    return lonlatalt;
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

  BBox2 GeoTransform::lonlat_to_lonlat_bbox(BBox2 const& lonlat_bbox) const {

    // Ensure we don't get incorrect results for empty boxes with strange corners.
    if (lonlat_bbox.empty())
      return BBox2();

    std::vector<vw::Vector2> points;
    sample_float_box(lonlat_bbox, points);

    BBox2 out_box;

    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try {
        out_box.grow(lonlat_to_lonlat(points[ptiter])); 
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
      }catch ( const std::exception & e ) {}
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
