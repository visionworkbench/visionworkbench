// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/Datum.h>
#include <vw/Math/Functions.h>

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

void vw::cartography::Datum::set_well_known_datum( std::string const& name ) {
  m_meridian_name = "Greenwich";
  m_geocentric = false;

  m_meridian_offset = 0;
  if (name == "WGS84" || name == "WGS_1984") {
    m_name = "WGS_1984";
    m_spheroid_name="WGS 84";
    m_semi_major_axis = 6378137.0;
    m_semi_minor_axis = 6356752.3;
    m_proj_str = "+ellps=WGS84 +datum=WGS84";
    return;
  }

  if (name == "WGS72" || name == "WGS_1972") {
    m_name="WGS_1972";
    m_spheroid_name="WGS 72";
    m_semi_major_axis = 6378135.0;
    m_semi_minor_axis = 6356750.5;
    m_proj_str = "+ellps=WGS72 +towgs84=0,0,4.5,0,0,0.554,0.2263";
    return;
  }

  if (name == "NAD83" || name == "North_American_Datum_1983") {
    m_name="North_American_Datum_1983";
    m_spheroid_name="GRS 1980";
    m_semi_major_axis = 6378137;
    m_semi_minor_axis = 6356752.3;
    m_proj_str = "+ellps=GRS80 +datum=NAD83";
    return;
  }

  if (name == "NAD27" || name == "North_American_Datum_1927") {
    m_name="North_American_Datum_1927";
    m_spheroid_name="Clarke 1866";
    m_semi_major_axis = 6378206.4;
    m_semi_minor_axis = 6356583.8;
    m_proj_str = "+ellps=clrk66 +datum=NAD27";
    return;
  }

  if (name == "D_MOON") {
    m_name = "D_MOON";
    m_spheroid_name = "MOON";
    m_meridian_name = "Reference Meridian";
    m_semi_major_axis = m_semi_minor_axis = 1737400;
    m_meridian_offset = 0.0;
    m_geocentric = false;
    m_proj_str = "+a=1737400 +b=1737400";
    return;
  }

  if (name == "D_MARS") {
    m_name = "D_MARS";
    m_spheroid_name = "MARS";
    m_meridian_name = "Reference Meridian";
    m_semi_major_axis = m_semi_minor_axis = 3396190;
    m_meridian_offset = 0.0;
    m_geocentric = false;
    m_proj_str = "+a=3396190 +b=3396190";
    return;
  }

  vw::vw_throw( vw::InputErr() << "Unknown datum string \"" << name << "\"!");
}

void vw::cartography::Datum::set_semi_major_axis(double val) {
  m_semi_major_axis = val;
  std::ostringstream strm;
  strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
  if (m_geocentric) strm << " +geoc";
  m_proj_str = strm.str();
}

void vw::cartography::Datum::set_semi_minor_axis(double val) {
  m_semi_minor_axis = val;
  std::ostringstream strm;
  strm << "+a=" << m_semi_major_axis << " +b=" << m_semi_minor_axis;
  if (m_geocentric) strm << " +geoc";
  m_proj_str = strm.str();
}

// return meridian radius of curvature.  NOT geocentric radius
double vw::cartography::Datum::radius(double /*lon*/, double lat) const {
  // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis;
  }

  // Bi-axial Ellpisoid datum
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

  // Bi-axial Ellpisoid datum
  // http://mathworld.wolfram.com/GeocentricLatitude.html
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
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
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;
  double slat = sin(M_PI/180*lat);
  return a / sqrt(1.0 - e2*slat*slat);
}

// return meridian radius of curvature.  NOT geocentric radius
double vw::cartography::Datum::geocentric_radius(double /*lon*/, double lat, double alt) const {
  // Optimize in the case of spherical datum
  if (m_semi_major_axis == m_semi_minor_axis) {
    return m_semi_major_axis + alt;
  }
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double Rn = a / sqrt(1.0-e2*slat*slat) + alt;

  return sqrt(Rn*Rn*(clat*clat + (1-e2)*(1-e2)*slat*slat));
}

double vw::cartography::Datum::inverse_flattening() const {
  return 1.0 / (1.0 - m_semi_minor_axis / m_semi_major_axis);
}

vw::Matrix3x3 vw::cartography::Datum::ecef_to_ned_matrix( vw::Vector3 const& p) const {
  double lat = p.y();
  if ( lat < -90 ) lat = -90;
  if ( lat > 90 ) lat = 90;

  double rlon = (p.x() + m_meridian_offset) * (M_PI/180);
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double slon = sin( rlon );
  double clon = cos( rlon );

  Matrix3x3 R;

  R(0,0) = -slat*clon;
  R(0,1) = -slat*slon;
  R(0,2) = clat;
  R(1,0) = -slon;
  R(1,1) = clon;
  R(1,2) = 0.0;
  R(2,0) = -clon*clat;
  R(2,1) = -slon*clat;
  R(2,2) = -slat;

  return R;
}


vw::Vector3 vw::cartography::Datum::geodetic_to_cartesian( vw::Vector3 const& p ) const {
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;

  double lat = p.y();
  if ( lat < -90 ) lat = -90;
  if ( lat > 90 ) lat = 90;

  double rlon = (p.x() + m_meridian_offset) * (M_PI/180);
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double slon = sin( rlon );
  double clon = cos( rlon );
  double radius = a / sqrt(1.0-e2*slat*slat);

  return Vector3( (radius+p.z()) * clat * clon,
                  (radius+p.z()) * clat * slon,
                  (radius*(1-e2)+p.z()) * slat );
}


// This function is based heavily on the similar
// Proj.4 function pj_Convert_Geocentric_To_Geodetic.
vw::Vector3 vw::cartography::Datum::cartesian_to_geodetic( vw::Vector3 const& p ) const {
  double a = m_semi_major_axis;
  double b = m_semi_minor_axis;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;

  static const double epsilon = 1.0e-12;
  static const double epsilon2 = epsilon*epsilon;
  static const int maxiter = 30;

  double normxy = sqrt(p.x()*p.x()+p.y()*p.y()); // distance between semi-minor axis and location
  double normp = norm_2(p);                      // distance between center and location

  double lon=0.0, alt=0.0;

  // compute the longitude
  if ( normxy/a < epsilon ) {
    // special case for the origin
    if ( normp/a < epsilon ) {
      return Vector3(0,90,-b);
    }
  }
  else {
    lon = atan2(p.y(),p.x()) / (M_PI/180);
  }

  // The following iterative algorithm was developped by
  // "Institut fur Erdmessung", University of Hannover, July 1988.
  // Internet: www.ife.uni-hannover.de
  double cgcl = normxy/normp;     // cos of geocentric latitude
  double sgcl = p.z()/normp;      // sin of geocentric latitude
  double rx = 1.0/sqrt(1.0-e2*(2.0-e2)*cgcl*cgcl);
  double clat = cgcl*(1.0-e2)*rx; // cos of geodetic latitude estimate
  double slat = sgcl*rx;          // sin of geodetic latitude estimate

  // loop to find lat (in quadrature) until |lat[i]-lat[i-1]|<epsilon, roughly
  for( int i=0; i<maxiter; ++i ) {
    double ri = a/sqrt(1.0-e2*slat*slat); // radius at estimated location
    alt = normxy*clat+p.z()*slat-ri*(1.0-e2*slat*slat);
    double rk = e2*ri/(ri+alt);
    rx = 1.0/sqrt(1.0-rk*(2.0-rk)*cgcl*cgcl);
    double new_clat = cgcl*(1.0-rk)*rx;
    double new_slat = sgcl*rx;
    // sin(lat[i]-lat[i-1]) ~= lat[i]-lat[i-1]
    double sdlat = new_slat*clat-new_clat*slat;
    clat = new_clat;
    slat = new_slat;
    if( sdlat*sdlat < epsilon2 ) break;
  }

  return Vector3( lon - m_meridian_offset, atan(slat/fabs(clat))/(M_PI/180), alt );
}

std::ostream& vw::cartography::operator<<( std::ostream& os, vw::cartography::Datum const& datum ) {
  os << "Geodeditic Datum --> Name: " << datum.name() << "  Spheroid: " << datum.spheroid_name()
     << "  Semi-major: " << datum.semi_major_axis()
     << "  Semi-minor: " << datum.semi_minor_axis()
     << "  Meridian: " << datum.meridian_name()
     << "  at " << datum.meridian_offset();
  return os;
}

