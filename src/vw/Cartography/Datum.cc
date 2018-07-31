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


#include <boost/algorithm/string.hpp>
#include <vw/Cartography/Datum.h>
#include <vw/Math/Functions.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>

vw::cartography::Datum::Datum(std::string const& name,
                              std::string const& spheroid_name,
                              std::string const& meridian_name,
                              double semi_major_axis,
                              double semi_minor_axis,
                              double meridian_offset)
  : m_name(name),
    m_spheroid_name(spheroid_name),
    m_meridian_name(meridian_name),
    m_semi_major_axis(semi_major_axis),
    m_semi_minor_axis(semi_minor_axis),
    m_meridian_offset(meridian_offset),
    m_geocentric(false)
{
  std::ostringstream strm;
  strm << "+a=" << semi_major_axis << " +b=" << semi_minor_axis;
  m_proj_str = strm.str();
}


// A wrapper around the GDAL/OGR API for setting the datum. Works for Earth datums.
void vw::cartography::Datum::set_datum_from_spatial_ref(OGRSpatialReference const& gdal_spatial_ref) {

  const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
  if (datum_name)
    this->name() = datum_name;

  const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
  if (spheroid_name)
    this->spheroid_name() = spheroid_name;

  const char* meridian_name = gdal_spatial_ref.GetAttrValue("PRIMEM");
  if (meridian_name)
    this->meridian_name() = meridian_name;

  OGRErr e1, e2;
  double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
  double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
  if (e1 != OGRERR_FAILURE && e2 != OGRERR_FAILURE) {
    this->set_semi_major_axis(semi_major);
    this->set_semi_minor_axis(semi_minor);
  }
  this->meridian_offset() = gdal_spatial_ref.GetPrimeMeridian();

  char* proj4_str_tmp;
  gdal_spatial_ref.exportToProj4(&proj4_str_tmp);

  this->proj4_str() = proj4_str_tmp;
  CPLFree( proj4_str_tmp );
}

// A wrapper around the GDAL/OGR API for setting the datum. Works for Earth datums.
void vw::cartography::Datum::set_datum_from_proj_str( std::string const& proj_str ) {

  OGRSpatialReference gdal_spatial_ref;
  if (gdal_spatial_ref.importFromProj4( proj_str.c_str() ))
    vw_throw( ArgumentErr() << "Failed to parse: \"" << proj_str << "\"." );

  set_datum_from_spatial_ref(gdal_spatial_ref);
  this->proj4_str() = proj_str; // The other call can change the string, don't let it!
}

void vw::cartography::Datum::set_well_known_datum( std::string const& name ) {
  m_meridian_name   = "Greenwich";
  m_geocentric      = false;
  m_meridian_offset = 0.0;

  // These numbers will be over-written later. However, we must
  // still initialize them, otherwise when the set_semi_major_axis()
  // function is invoked later it will result in un-initialized
  // variables (since that function sets the semi-major axis but
  // assumes the semi-minor axis is already set).
  m_semi_major_axis =  6378137;
  m_semi_minor_axis = 6356752.3142;

  std::string up_name = boost::to_upper_copy(name);

  if (up_name == "WGS84"    || up_name == "WGS_1984" ||
      up_name == "WGS 1984" || up_name == "WGS1984"   ||
      up_name == "WORLD GEODETIC SYSTEM 1984" || up_name == "EARTH") {
    set_datum_from_proj_str("+proj=longlat +datum=WGS84 +no_defs");
    return;
  }

  if (up_name == "WGS72" || up_name == "WGS_1972") {
    set_datum_from_proj_str("+proj=longlat +ellps=WGS72 +no_defs");
    return;
  }

  if (up_name == "NAD83" ||
      up_name == boost::to_upper_copy(std::string("North_American_Datum_1983"))) {
    set_datum_from_proj_str("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs");
    return;
  }

  if (up_name == "NAD27" ||
      up_name == boost::to_upper_copy(std::string("North_American_Datum_1927"))) {
    set_datum_from_proj_str("+proj=longlat +datum=NAD27 +no_defs");
    return;
  }

  if (up_name == "D_MOON" || up_name == "MOON") {
    m_name            = "D_MOON";
    m_spheroid_name   = "MOON";
    m_meridian_name   = "Reference Meridian";
    m_semi_major_axis = m_semi_minor_axis = 1737400;
    m_meridian_offset = 0.0;
    m_proj_str        = "+a=1737400 +b=1737400";
    return;
  }

  if (up_name == "D_MARS" || up_name == "MARS") {
    m_name            = "D_MARS";
    m_spheroid_name   = "MARS";
    m_meridian_name   = "Reference Meridian";
    m_semi_major_axis = m_semi_minor_axis = 3396190;
    m_meridian_offset = 0.0;
    m_proj_str        = "+a=3396190 +b=3396190";
    return;
  }

  if (up_name == "MOLA") {
    m_name            = "D_MARS";
    m_spheroid_name   = "MARS";
    m_meridian_name   = "Reference Meridian";
    m_semi_major_axis = m_semi_minor_axis = 3396000;
    m_meridian_offset = 0.0;
    m_proj_str        = "+a=3396000 +b=3396000";
    return;
  }

  vw::vw_throw( vw::InputErr() << "Unknown datum: " << name << ".");
}

void vw::cartography::Datum::set_semi_major_axis(double val) {
  m_semi_major_axis = val;
  std::ostringstream strm;
  strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
  if (m_geocentric)
    strm << " +geoc";
  m_proj_str = strm.str();
}

void vw::cartography::Datum::set_semi_minor_axis(double val) {
  m_semi_minor_axis = val;
  std::ostringstream strm;
  strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
  if (m_geocentric)
    strm << " +geoc";
  m_proj_str = strm.str();
}

// return meridian radius of curvature.  NOT geocentric radius
double vw::cartography::Datum::radius(double /*lon*/, double lat) const {
  // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis;
  }

  // Bi-axial ellpisoid datum
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
  double t = atan((a/b) * tan(lat * M_PI / 180.0));
  double x = a * cos(t);
  double y = b * sin(t);
  return sqrt(x*x + y*y);
}

double vw::cartography::Datum::geocentric_latitude(double lat) const {
   // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis;
  }

  // Bi-axial ellpisoid datum
  // http://mathworld.wolfram.com/GeocentricLatitude.html
  double a  = m_semi_major_axis;
  double b  = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;
  return atan((1-e2)*tan(lat * M_PI / 180.0));

}

double vw::cartography::Datum::radius_of_curvature(double /*lon*/, double lat) const {
  // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis;
  }

  // Bi-axial Ellpisoid datum
  double a    = m_semi_major_axis;
  double b    = m_semi_minor_axis;
  double a2   = a * a;
  double b2   = b * b;
  double e2   = (a2 - b2) / a2;
  double slat = sin(M_PI/180*lat);
  return a / sqrt(1.0 - e2*slat*slat);
}

// return meridian radius of curvature.  NOT geocentric radius
double vw::cartography::Datum::geocentric_radius(double /*lon*/, double lat, double alt) const {
  // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis + alt;
  }
  double a    = m_semi_major_axis;
  double b    = m_semi_minor_axis;
  double a2   = a * a;
  double b2   = b * b;
  double e2   = (a2 - b2) / a2;
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double Rn   = a / sqrt(1.0-e2*slat*slat) + alt;

  return sqrt(Rn*Rn*(clat*clat + (1-e2)*(1-e2)*slat*slat));
}

double vw::cartography::Datum::inverse_flattening() const {
  return 1.0 / (1.0 - m_semi_minor_axis / m_semi_major_axis);
}

vw::Matrix3x3 vw::cartography::Datum::lonlat_to_ned_matrix( vw::Vector2 const& lonlat) const {
  double lon = lonlat.x();
  double lat = lonlat.y();
  if ( lat < -90 ) lat = -90;
  if ( lat >  90 ) lat =  90;

  double rlon = (lon + m_meridian_offset) * (M_PI/180);
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double slon = sin( rlon );
  double clon = cos( rlon );

  Matrix3x3 R;

  R(0,0) = -slat*clon;
  R(1,0) = -slat*slon;
  R(2,0) = clat;
  R(0,1) = -slon;
  R(1,1) = clon;
  R(2,1) = 0.0;
  R(0,2) = -clon*clat;
  R(1,2) = -slon*clat;
  R(2,2) = -slat;

  return R;
}

vw::Vector3 vw::cartography::Datum::geodetic_to_cartesian( vw::Vector3 const& llh ) const {
  double a  = m_semi_major_axis;
  double b  = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;

  double lat = llh.y();
  if ( lat < -90 ) lat = -90;
  if ( lat >  90 ) lat = 90;

  double rlon = (llh.x() + m_meridian_offset) * (M_PI/180);
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double slon = sin( rlon );
  double clon = cos( rlon );
  double radius = a / sqrt(1.0-e2*slat*slat);

  return Vector3( (radius+llh.z()) * clat * clon,
                  (radius+llh.z()) * clat * slon,
                  (radius*(1-e2)+llh.z()) * slat );
}

// This algorithm is a non-iterative algorithm from "An analytical
// method to transform geocentric into geodetic coordinates" by Hugues
// Vermeille, Journal of Geodesy 2011.
//
// This is an improvement over the 1988/Proj4's implementation as it's
// a smidgen faster and it still works near the center of the datum.
vw::Vector3 vw::cartography::Datum::cartesian_to_geodetic( vw::Vector3 const& xyz ) const {
  const double a2 = m_semi_major_axis * m_semi_major_axis;
  const double b2 = m_semi_minor_axis * m_semi_minor_axis;
  const double e2 = 1 - b2 / a2;
  const double e4 = e2 * e2;

  double xy_dist = sqrt( xyz[0] * xyz[0] + xyz[1] * xyz[1] );
  double p  = ( xyz[0] * xyz[0] + xyz[1] * xyz[1] ) / a2;
  double q  = ( 1 - e2 ) * xyz[2] * xyz[2] / a2;
  double r  = ( p + q - e4 ) / 6.0;
  double r3 = r * r * r;

  Vector3 llh;

  double evolute = 8 * r3 + e4 * p * q;
  double u = std::numeric_limits<double>::quiet_NaN();
  if ( evolute > 0 ) {
    // outside the evolute
    double right_inside_pow = sqrt(e4 * p * q);
    double sqrt_evolute = sqrt( evolute );
    u = r + 0.5 * pow(sqrt_evolute + right_inside_pow, 2.0/3.0)
          + 0.5 * pow(sqrt_evolute - right_inside_pow, 2.0/3.0);
  } else if ( fabs(xyz[2]) < std::numeric_limits<double>::epsilon() ) {
    // On the equator plane
    llh[1] = 0;
    llh[2] = norm_2( xyz ) - m_semi_major_axis;
  } else if ( evolute < 0 and fabs(q) > std::numeric_limits<double>::epsilon() ) {
    // On or inside the evolute
    double atan_result = atan2( sqrt( e4 * p * q ), sqrt( -evolute ) + sqrt(-8 * r3) );
    u = -4 * r * sin( 2.0 / 3.0 * atan_result )
               * cos( M_PI / 6.0 + 2.0 / 3.0 * atan_result );
  } else if ( fabs(q) < std::numeric_limits<double>::epsilon() and p <= e4 ) {
    // In the singular disc
    llh[2] = -m_semi_major_axis * sqrt(1 - e2) * sqrt(e2 - p) / sqrt(e2);
    llh[1] = 2 * atan2( sqrt(e4 - p), sqrt(e2*(e2 - p)) + sqrt(1-e2) * sqrt(p) );
  } else {
    // Near the cusps of the evolute
    double inside_pow = sqrt(evolute) + sqrt(e4 * p * q);
    u = r + 0.5   * pow(inside_pow, 2.0/3.0)
          + 2*r*r * pow(inside_pow,-2.0/3.0);
  }

  if (!std::isnan(u) ) {
    double v   = sqrt( u * u + e4 * q );
    double u_v = u + v;
    double w   = e2 * ( u_v - q ) / ( 2 * v );
    double k   = u_v / ( w + sqrt( w * w + u_v ) );
    double D   = k * xy_dist / ( k + e2 );
    double dist_2 = D * D + xyz[2] * xyz[2];
    llh[2] = ( k + e2 - 1 ) * sqrt( dist_2 ) / k;
    llh[1] = 2 * atan2( xyz[2], sqrt( dist_2 ) + D );
  }

  if ( xy_dist + xyz[0] > ( sqrt(2) - 1 ) * xyz[1] ) {
    // Longitude is between -135 and 135
    llh[0] = 360.0 * atan2( xyz[1], xy_dist + xyz[0] ) / M_PI;
  } else if ( xy_dist + xyz[1] < ( sqrt(2) + 1 ) * xyz[0] ) {
    // Longitude is between -225 and 45
    llh[0] = - 90.0 + 360.0 * atan2( xyz[0], xy_dist - xyz[1] ) / M_PI;
  } else {
    // Longitude is between -45 and 225
    llh[0] = 90.0 - 360.0 * atan2( xyz[0], xy_dist + xyz[1] ) / M_PI;
  }
  llh[0] -= m_meridian_offset;
  llh[1] *= 180.0 / M_PI;

  return llh;
}

std::ostream& vw::cartography::operator<<( std::ostream& os, vw::cartography::Datum const& datum ) {
  std::ostringstream oss; // To use custom precision
  oss.precision(17);
  oss << "Geodetic Datum --> Name: " << datum.name() << "  Spheroid: " << datum.spheroid_name()
      << "  Semi-major axis: " << datum.semi_major_axis()
      << "  Semi-minor axis: " << datum.semi_minor_axis()
      << "  Meridian: "   << datum.meridian_name() << " at " << datum.meridian_offset()
      << "  Proj4 Str: "  << datum.proj4_str();
  os << oss.str();
  return os;
}

// Free associated functions

vw::Vector3
vw::cartography::datum_intersection(double semi_major_axis, double semi_minor_axis,
                                    vw::Vector3 camera_ctr, vw::Vector3 camera_vec) {

  // The datum is a spheroid. To simplify the calculations, scale
  // everything in such a way that the spheroid becomes a
  // sphere. Scale back at the end of computation.

  double z_scale = semi_major_axis / semi_minor_axis;
  camera_ctr.z() *= z_scale;
  camera_vec.z() *= z_scale;
  camera_vec = normalize(camera_vec);
  double radius_2 = semi_major_axis * semi_major_axis;
  double alpha = -dot_prod(camera_ctr, camera_vec );
  vw::Vector3 projection = camera_ctr + alpha*camera_vec;
  if ( norm_2_sqr(projection) > radius_2 ) {
    // did not intersect
    return vw::Vector3();
  }

  alpha -= sqrt( radius_2 -
                 norm_2_sqr(projection) );
  vw::Vector3 intersection = camera_ctr + alpha * camera_vec;
  intersection.z() /= z_scale;
  return intersection;
}

// Intersect the ray back-projected from the camera with the datum.
vw::Vector3
vw::cartography::datum_intersection( vw::cartography::Datum const& datum,
                                     vw::Vector3 camera_ctr, vw::Vector3 camera_vec) {
  return vw::cartography::datum_intersection(datum.semi_major_axis(), datum.semi_minor_axis(),
                                             camera_ctr, camera_vec);
}
