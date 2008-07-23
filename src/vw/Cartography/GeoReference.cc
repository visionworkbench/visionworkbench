// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#include <vw/Cartography/GeoReference.h>

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#endif

#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

// Boost
#include <boost/algorithm/string.hpp>

// Proj.4
#include <projects.h>


void vw::cartography::read_georeference( vw::cartography::GeoReference& georef, vw::ImageResource const& resource ) {
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  DiskImageResourceGDAL const* gdal = dynamic_cast<DiskImageResourceGDAL const*>( &resource );
  if( gdal ) return read_gdal_georeference( georef, *gdal );
#endif
  DiskImageResourcePDS const* pds = dynamic_cast<DiskImageResourcePDS const*>( &resource );
  if( pds ) return read_pds_georeference( georef, *pds );
  vw_throw(NoImplErr() << "This image resource does not support reading georeferencing information.");
}

void vw::cartography::write_georeference( vw::ImageResource& resource, vw::cartography::GeoReference const& georef ) {
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  DiskImageResourceGDAL* gdal = dynamic_cast<DiskImageResourceGDAL*>( &resource );
  if( gdal ) return write_gdal_georeference( *gdal, georef );
#endif
  // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
  vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
}

namespace vw {
namespace cartography {


  std::string GeoReference::proj4_str() const { 
    return m_proj_projection_str;
  }

  std::string GeoReference::overall_proj4_str() const {
    std::string proj4_str = m_proj_projection_str + " " + m_datum.proj4_str() + " +no_defs";
    return proj4_str;
  }
  
  void GeoReference::init_proj() {
    m_proj_context = boost::shared_ptr<ProjContext>(new ProjContext(overall_proj4_str()));
  }

  /// Construct a default georeference.  This georeference will use
  /// the identity matrix as the initial transformation matrix, and
  /// select the default datum (WGS84) and projection (geographic).
  GeoReference::GeoReference() {
    set_transform(vw::math::identity_matrix<3>());
    set_geographic();
    init_proj();
  }

  /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(Datum const& datum) : GeoReferenceBase(datum){
    set_transform(vw::math::identity_matrix<3>());
    set_geographic();
    init_proj();
  }
  
  /// Takes a geodetic datum and an affine transformation matrix
  GeoReference::GeoReference(Datum const& datum, Matrix<double,3,3> const& transform) : GeoReferenceBase(datum) {
    set_transform(transform);
    set_geographic();
    init_proj();
  }

  void GeoReference::set_transform(Matrix3x3 transform) {
    m_transform = transform;
    m_shifted_transform = m_transform;
    m_shifted_transform(0,2) += 0.5*m_transform(0,0);
    m_shifted_transform(1,2) += 0.5*m_transform(1,1);
    m_inv_transform = vw::math::inverse(m_transform);
    m_inv_shifted_transform = vw::math::inverse(m_shifted_transform);
  }

  // We override the base classes method here so that we have the
  // opportunity to call init_proj()
  void GeoReference::set_datum(Datum const& datum) {
    m_datum = datum;
    init_proj();
  }
  
  // Adjust the affine transform to the VW convention ( [0,0] is at
  // the center of upper left pixel) if file is georeferenced
  // according to the convention that [0,0] is the upper left hand
  // corner of the upper left pixel.
  inline Matrix3x3 const& GeoReference::vw_native_transform() const {
    if (m_pixel_interpretation == GeoReference::PixelAsArea)
      return m_shifted_transform;
    else
      return m_transform;
  }

  inline Matrix3x3 const& GeoReference::vw_native_inverse_transform() const {
    if (m_pixel_interpretation == GeoReference::PixelAsArea) 
      return m_inv_shifted_transform;
    else
      return m_inv_transform;
  }

  void GeoReference::set_well_known_geogcs(std::string name) {
    m_datum.set_well_known_datum(name);
    init_proj();
  }

  void GeoReference::set_geographic() {
    m_is_projected = false;
    m_proj_projection_str = "+proj=longlat";
    init_proj();
  }

  void GeoReference::set_sinusoidal(double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sinu +lon_0=" << center_longitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }
  
  void GeoReference::set_mercator(double center_latitude, double center_longitude, double latitude_of_true_scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=merc +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +lat_ts=" << latitude_of_true_scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=tmerc +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_orthographic(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=ortho +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_oblique_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sterea +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=stere +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=laea +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_lambert_conformal(double std_parallel_1, double std_parallel_2, double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=lcc +lat_1=" << std_parallel_1 << " +lat_2=" << std_parallel_2 << " +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }
  
  void GeoReference::set_UTM(int zone, int north) {
    std::ostringstream strm;
    strm << "+proj=utm +zone=" << zone;
    if (!north) strm << " +south";
    strm << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_proj4_projection_str(std::string const& s) { 
    m_proj_projection_str = s; 
    if (s.find("+proj=longlat") == 0)
      m_is_projected = false;
    else
      m_is_projected = true;
    init_proj();
  }

  /// For a given pixel coordinate, compute the position of that
  /// pixel in this georeferenced space.
  Vector2 GeoReference::pixel_to_point(Vector2 pix) const {
    Vector2 loc;
    Matrix<double,3,3> M = this->vw_native_transform();
    double denom = pix[0] * M(2,0) + pix[1] * M(2,1) + M(2,2);
    loc[0] = (pix[0] * M(0,0) + pix[1] * M(0,1) + M(0,2)) / denom;
    loc[1] = (pix[0] * M(1,0) + pix[1] * M(1,1) + M(1,2)) / denom;
    return loc;
  }
  
  /// For a given location 'loc' in projected space, compute the
  /// corresponding pixel coordinates in the image.
  Vector2 GeoReference::point_to_pixel(Vector2 loc) const {
    Vector2 pix;
    Matrix<double,3,3> M = this->vw_native_inverse_transform();
    double denom = loc[0] * M(2,0) + loc[1] * M(2,1) + M(2,2);
    pix[0] = (loc[0] * M(0,0) + loc[1] * M(0,1) + M(0,2)) / denom;
    pix[1] = (loc[0] * M(1,0) + loc[1] * M(1,1) + M(1,2)) / denom;
    return pix;
  }

  
  /// For a point in the projected space, compute the position of
  /// that point in unprojected (Geographic) coordinates (lat,lon).
  Vector2 GeoReference::point_to_lonlat(Vector2 loc) const {
    if ( ! m_is_projected ) return loc;

    XY projected;  
    LP unprojected;

    projected.u = loc[0];
    projected.v = loc[1];

    unprojected = pj_inv(projected, m_proj_context->proj_ptr());
    CHECK_PROJ_ERROR;

    // Convert from radians to degrees.
    return Vector2(unprojected.u * RAD_TO_DEG, unprojected.v * RAD_TO_DEG);
  }
  
  /// Given a position in geographic coordinates (lat,lon), compute
  /// the location in the projected coordinate system.
  Vector2 GeoReference::lonlat_to_point(Vector2 lon_lat) const {
    if ( ! m_is_projected ) return lon_lat;

    // Clamp the latitude range to [-90, 90] as occasionally we get edge 
    // pixels that extend slightly beyond that range and cause Proj.4 to 
    // fail.
    if(lon_lat[1] > 90) lon_lat[1] = 90;
    else if(lon_lat[1] < -90) lon_lat[1] = -90;

    XY projected;  
    LP unprojected;

    // Proj.4 expects the (lon,lat) pair to be in radians, so we
    // must make a conversion if the CS in geographic (lat/lon).
    unprojected.u = lon_lat[0] * DEG_TO_RAD;
    unprojected.v = lon_lat[1] * DEG_TO_RAD;

    projected = pj_fwd(unprojected, m_proj_context->proj_ptr());
    // Needed because latitudes -90 and 90 can fail on spheroids.
    if(pj_errno == -20)
      if(lon_lat[1] == 90 || lon_lat[1] == -90)
        return Vector2(-1, -1);
    CHECK_PROJ_ERROR;

    return Vector2(projected.u, projected.v);
  }

  /************** Functions for class ProjContext *******************/
  char** ProjContext::split_proj4_string(std::string const& proj4_str, int &num_strings) {
    std::vector<std::string> arg_strings;
    std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
    boost::split( arg_strings, trimmed_proj4_str, boost::is_any_of(" ") ); 

    char** strings = new char*[arg_strings.size()]; 
    for ( unsigned i = 0; i < arg_strings.size(); ++i ) {
      strings[i] = new char[2048];
      strncpy(strings[i], arg_strings[i].c_str(), 2048);
    }
    num_strings = arg_strings.size();
    return strings;
  }

  ProjContext::ProjContext(std::string const& proj4_str) {

    // proj.4 is expecting the parameters to be split up into seperate
    // c-style strings.
    int num;
    char** proj_strings = split_proj4_string(proj4_str, num);
    m_proj_ptr = pj_init(num, proj_strings);
    CHECK_PROJ_INIT_ERROR(proj4_str);

    for (int i = 0; i < num; i++) 
      delete [] proj_strings[i];
    delete [] proj_strings;
  }

  ProjContext::~ProjContext() {
    pj_free(m_proj_ptr);
  }

  /***************** Functions for output GeoReferences *****************/
  namespace output {  
    GeoReference kml::get_output_georeference(int xresolution, int yresolution) {
      GeoReference r;
      Matrix3x3 transform;

      r.set_well_known_geogcs("WGS84");

      // We specify here the KML transformation for longitude/latitude
      // to an exact pixel. Thus when we set the georeference 
      // transform (which is pixel -> lon/lat), we have to use its 
      // inverse.
      transform(0,0) = xresolution / 360.0;
      transform(0,1) = 0;
      transform(0,2) = xresolution / 2.0 - 0.5;
      transform(1,0) = 0;
      transform(1,1) = -(yresolution / 360.0);
      transform(1,2) = yresolution / 2.0 - 0.5;
      transform(2,0) = 0;
      transform(2,1) = 0;
      transform(2,2) = 1;
      r.set_transform(vw::math::inverse(transform));

      return r;
    }

    GeoReference tms::get_output_georeference(int resolution) {
      GeoReference r;
      Matrix3x3 transform;

      r.set_well_known_geogcs("WGS84");
      // We specify here the TMS transformation for longitude/latitude
      // to an exact pixel. Thus when we set the georeference 
      // transform (which is pixel -> lon/lat), we have to use its 
      // inverse.
      transform(0,0) = resolution / 360.0;
      transform(0,1) = 0;
      transform(0,2) = resolution / 2.0 - 0.5;
      transform(1,0) = 0;
      transform(1,1) = -(resolution / 360.0);
      transform(1,2) = .75 * resolution - 0.5;
      transform(2,0) = 0;
      transform(2,1) = 0;
      transform(2,2) = 1;
      r.set_transform(vw::math::inverse(transform));

      return r;
    }

  } // namespace vw::cartography::output
}} // vw::cartography
