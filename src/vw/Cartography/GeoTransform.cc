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

// TODO(oalexan1): Wipe proj
#include <proj.h>

namespace vw {
namespace cartography {

  // Wipe the +proj and other entries from a projection before appending a new one,
  // as otherwise PROJ gets confused. Must make sure to keep the radius.
  // TODO(oalexan1): This is a hack. Need to figure out how to properly
  // do conversions between datums. 
  // TODO(oalexan1): Wipe this
  std::string wipe_proj(std::string const& str_in) {

    std::vector<std::string> tokens;
    boost::split(tokens, boost::trim_copy(str_in), boost::is_any_of(" "));

    std::vector<std::string> to_wipe = {"+proj", "+lon", "+lat", "+k", "+x", "+y"};
    
    std::string str_out = "";
    for (size_t k = 0; k < tokens.size(); k++) {
      
      // Wipe any tokens that start with the strings in to_wipe
      bool will_skip = false;
      for (size_t i = 0; i < to_wipe.size(); i++) {
        if (boost::starts_with(tokens[k], to_wipe[i])) {
          will_skip = true;
          break;
        }
      }
      if (will_skip)
        continue;

      if (!str_out.empty())
        str_out += " ";

      str_out += tokens[k];
    }

    return str_out;
  }
  
  using vw::math::BresenhamLine;

  // Create transform between given proj strings. Note that what is returned
  // are pointers, by reference.
  void create_proj_transform(std::string const& src_proj_str,
                             std::string const& dst_proj_str,
                             PJ_CONTEXT* & pj_context,
                             PJ*         & pj_transform) {

    pj_context = proj_context_create();
    PJ_AREA *area = NULL;
    pj_transform = proj_create_crs_to_crs(pj_context, src_proj_str.c_str(),
                                          dst_proj_str.c_str(), area);
    if (pj_transform == NULL)
      vw::vw_throw(vw::ArgumentErr() 
                   << "Failed to initialize a geotransform for PROJ strings: "
                   << src_proj_str << " and " 
                   << dst_proj_str << ".\n");
  }
  
  // Constructor
  GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                             BBox2 const& src_bbox, BBox2 const& dst_bbox) :
    m_src_georef(src_georef), m_dst_georef(dst_georef),
    m_src_bbox(src_bbox), m_dst_bbox(dst_bbox) {
     
    // This are NULL for now, may be initialized if needed to convert among datums
    m_pj_context = NULL;
    m_pj_transform = NULL;
 
    const std::string src_datum = m_src_georef.datum().get_wkt();
    const std::string dst_datum = m_dst_georef.datum().get_wkt();

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
    
    // This optimizes the case where the two datums are the same,
    // and thus we don't need to call proj to convert between them
    // as we transform.
    if (src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {
      m_skip_datum_conversion = false;

      // Set up longlat coordinate systems with given datums
      
      // TODO(oalexan1): Set here honest CRS that use the lonlat projection
      // wipe the PROJ.4 logic
      
      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the
      // source or destination georef uses.
      ss_src << "+proj=longlat " << wipe_proj(src_datum);
      m_src_datum_proj_str = ss_src.str();
      
      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=longlat " << wipe_proj(dst_datum);
      m_dst_datum_proj_str = ss_dst.str();
      
      // This transform will only be used to convert lon,lat,alt between
      // datums. It will not be used to convert between projections.
      create_proj_transform(m_src_datum_proj_str, m_dst_datum_proj_str,  
                            m_pj_context, m_pj_transform); // outputs
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
      double alt = 0.0;
      double max_err = 0.0;
      for (double lon = -180; lon <= 180; lon += 90) {
        for (double lat = -90; lat <= 90; lat += 45) {
          Vector2 lon_lat(lon, lat);
          bool forward = true;
          Vector2 lon_lat2 = GeoTransform::lonlat_to_lonlat(lon_lat, forward);
          max_err = std::max(max_err, norm_2(lon_lat - lon_lat2));
          Vector2 lon_lat3 = GeoTransform::lonlat_to_lonlat(lon_lat2, !forward);
          max_err = std::max(max_err, norm_2(lon_lat2 - lon_lat3));
        }
      }
      
      if (max_err < 1e-8)
        m_skip_datum_conversion = true;
    }
    
    // Cannot reliably transform between images with different datums,
    // as that requires a vertical datum shift, while for GeoTransform
    // we assume a purely 2D transform.
    if (!m_skip_datum_conversion)
      vw::vw_throw(vw::ArgumentErr() << "Cannot handle images with different datums. "
                   << "Use gdalwarp to convert them to the same datum.\n");
  }

  GeoTransform::GeoTransform(GeoTransform const& other) {
    *this = other;
  }

  GeoTransform& GeoTransform::operator=(GeoTransform const& other) {
    m_src_georef            = other.m_src_georef;
    m_dst_georef            = other.m_dst_georef;
    m_src_bbox              = other.m_src_bbox;
    m_dst_bbox              = other.m_dst_bbox;
    m_src_datum_proj_str    = other.m_src_datum_proj_str;
    m_dst_datum_proj_str    = other.m_dst_datum_proj_str;
    m_skip_map_projection   = other.m_skip_map_projection;
    m_skip_datum_conversion = other.m_skip_datum_conversion;

    // A mutex seems necessary to avoid a crash. Presumably
    // this creation logic does not like to be created from multiple threads.
    // TODO(oalexan1): It is not clear if this is either necessary
    // or sufficient to avoid the crash.
    if (!m_skip_datum_conversion) {
      Mutex::WriteLock write_lock(m_mutex);
      create_proj_transform(m_src_datum_proj_str, m_dst_datum_proj_str,  
                            m_pj_context, m_pj_transform); // outputs
    }
    
    return *this;
  }

  GeoTransform::~GeoTransform() {
    // TODO(oalexan1): Must deallocate pj_transform and pj_context. But for now
    // don't have a testcase, as GeoTransform does not initialize these given
    // that datum conversion is not allowed.
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
    bool forward = false;
    if (m_skip_map_projection)
      return m_src_georef.point_to_pixel(m_dst_georef.pixel_to_point(v));
    Vector2 dst_lonlat = m_dst_georef.pixel_to_lonlat(v);
    
    // For a georef with projection, must ensure input lonlat is in the
    // [-180,180] range.
    if (m_src_georef.is_projected())
      wrapLon(dst_lonlat[0]);
    
    if (m_skip_datum_conversion)
      return m_src_georef.lonlat_to_pixel(dst_lonlat);
     
     // This will require a datum conversion, which is slow
     Vector2 src_lonlat = lonlat_to_lonlat(dst_lonlat, forward);
     return m_src_georef.lonlat_to_pixel(src_lonlat);
  }

  /// Given a pixel coordinate of an image in a source
  /// georeference frame, this routine computes the corresponding
  /// pixel in the destination (transformed) image.
  Vector2 GeoTransform::forward(Vector2 const& v) const {
    bool forward = true;
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(m_src_georef.pixel_to_point(v));

    Vector2 src_lonlat = m_src_georef.pixel_to_lonlat(v);
    
    // For a geo ref with projection, must ensure input lonlat is in the
    // [-180,180] range.
    if (m_dst_georef.is_projected())
      wrapLon(src_lonlat[0]);
    
    if (m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_pixel(src_lonlat);
    
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, forward);
    return m_dst_georef.lonlat_to_pixel(dst_lonlat);
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
    if (m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_point(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_point(dst_lonlat);
  }


  Vector2 GeoTransform::pixel_to_point(Vector2 const& v) const {
    Vector2 src_point = m_src_georef.pixel_to_point(v);
    if (m_skip_map_projection)
      return src_point;
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(src_point);
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_point(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_point(dst_lonlat);
  }

  Vector2 GeoTransform::point_to_pixel(Vector2 const& v) const {
    if (m_skip_map_projection)
      return m_dst_georef.point_to_pixel(v);
     
    Vector2 src_lonlat = m_src_georef.point_to_lonlat(v);
    if (m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_pixel(src_lonlat);
    Vector2 dst_lonlat = lonlat_to_lonlat(src_lonlat, true);
    return m_dst_georef.lonlat_to_pixel(dst_lonlat);
  }

  // Performs a forward or reverse datum conversion.
  // TODO(oalexan1): Must consider 360 degree offsets here.
  // TODO(oalexan1): Must consider altitude here always!
  Vector2 GeoTransform::lonlat_to_lonlat(Vector2 const& lonlat, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlat;

    // Note that the vertical component of the transform is not used
    Vector3 lonlatalt(lonlat[0], lonlat[1], 0.0);
    return subvector(GeoTransform::lonlatalt_to_lonlatalt(lonlatalt, forward), 0, 2);
  }

  // Note: Not sure if this is ever reached as the datum conversion is not supported.
  Vector3 GeoTransform::lonlatalt_to_lonlatalt(Vector3 const& lonlatalt, bool forward) const {
    if (m_skip_datum_conversion)
      return lonlatalt;

    double lon_rad = proj_torad(lonlatalt[0]); // proj4 requires radians
    double lat_rad = proj_torad(lonlatalt[1]);
    double alt = lonlatalt[2];

    PJ_COORD c_in, c_out;
    c_in.lpzt.t = HUGE_VAL;
    c_in.lpzt.z = alt;
    c_in.lp.lam = lon_rad;
    c_in.lp.phi = lat_rad;

    // TODO(oalexan1): Not clear about the mutex
    Mutex::WriteLock write_lock(m_mutex);
    if (forward) 
      c_out = proj_trans(m_pj_transform, PJ_FWD, c_in);
    else
      c_out = proj_trans(m_pj_transform, PJ_INV, c_in);
    
    vw::Vector3 out(proj_todeg(c_out.lp.lam), proj_todeg(c_out.lp.phi), c_out.lpzt.z);
    return out;
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

