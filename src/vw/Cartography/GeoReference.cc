// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Cartography/GeoReference.h>

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include "ogr_spatialref.h"
#include "cpl_string.h"
#endif

#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

// Boost
#include <boost/algorithm/string.hpp>

// Proj.4
#include <projects.h>


bool vw::cartography::read_georeference( vw::cartography::GeoReference& georef,
                                         vw::ImageResource const& resource ) {

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  DiskImageResourceGDAL const* gdal =
    dynamic_cast<DiskImageResourceGDAL const*>( &resource );
  if( gdal ) return read_gdal_georeference( georef, *gdal );
#endif

  DiskImageResourcePDS const* pds =
    dynamic_cast<DiskImageResourcePDS const*>( &resource );
  if( pds ) return read_pds_georeference( georef, *pds );
  return false;
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

  void GeoReference::set_equirectangular(double center_latitude, double center_longitude, double latitude_of_true_scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=eqc +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +lat_ts=" << latitude_of_true_scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
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

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL
  void GeoReference::set_wkt(std::string const& wkt) {
    const char *wkt_str = wkt.c_str();
    char **wkt_ptr = (char**)(&wkt_str);

    OGRSpatialReference gdal_spatial_ref;
    gdal_spatial_ref.importFromWkt(wkt_ptr);

    // Read projection information out of the file
    char* proj_str_tmp;
    gdal_spatial_ref.exportToProj4(&proj_str_tmp);
    std::string proj4_str = proj_str_tmp;
    CPLFree( proj_str_tmp );
    // For debugging:
    //      vw_out() << "PROJ in --> " << proj4_str << "\n";

    std::vector<std::string> input_strings;
    std::vector<std::string> output_strings;
    std::vector<std::string> datum_strings;
    std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
    boost::split( input_strings, trimmed_proj4_str, boost::is_any_of(" ") );
    for (unsigned int i = 0; i < input_strings.size(); ++i) {

      // Pick out the parts of the projection string that pertain to
      // map projections.  We essentially want to eliminate all of
      // the strings that have to do with the datum, since those are
      // handled by interacting directly with the
      // OGRSpatialReference below. This is sort of messy, but it's
      // the easiest way to do this, as far as I can tell.
      if ((input_strings[i].find("+proj=") == 0) ||
          (input_strings[i].find("+x_0=") == 0) ||
          (input_strings[i].find("+y_0=") == 0) ||
          (input_strings[i].find("+lon") == 0) ||
          (input_strings[i].find("+lat") == 0) ||
          (input_strings[i].find("+k=") == 0) ||
          (input_strings[i].find("+lat_ts=") == 0) ||
          (input_strings[i].find("+ns") == 0) ||
          (input_strings[i].find("+no_cut") == 0) ||
          (input_strings[i].find("+h=") == 0) ||
          (input_strings[i].find("+W=") == 0) ||
          (input_strings[i].find("+units=") == 0) ||
          (input_strings[i].find("+zone=") == 0)) {
        output_strings.push_back(input_strings[i]);
      } else if ((input_strings[i].find("+ellps=") == 0) ||
                 (input_strings[i].find("+datum=") == 0)) {
        // We put these in the proj4_str for the Datum class.
        datum_strings.push_back(input_strings[i]);
      }
    }
    std::ostringstream strm;
    for (unsigned int i = 0; i < output_strings.size(); ++i) {
      strm << output_strings[i] << " ";
    }
    // For debugging:
    //      vw_out() << "     out --> " << strm.str() << "\n";

    // If the file contains no projection related information, we
    // supply proj.4 with a "default" interpretation that the file
    // is in geographic (unprojected) coordinates.
    if (output_strings.empty())
      set_proj4_projection_str("+proj=longlat");
    else
      set_proj4_projection_str(strm.str());

    int utm_north = 0;
    int utm_zone = gdal_spatial_ref.GetUTMZone(&utm_north);
    if (utm_zone)
      set_UTM(utm_zone, utm_north);

    // Read in the datum information
    Datum datum;
    const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
    if (datum_name) { datum.name() = datum_name; }
    const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
    if (spheroid_name) { datum.spheroid_name() = spheroid_name; }
    const char* meridian_name = gdal_spatial_ref.GetAttrValue("PRIMEM");
    if (meridian_name) { datum.meridian_name() = meridian_name; }
    OGRErr e1, e2;
    double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
    double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
    if (e1 != OGRERR_FAILURE && e2 != OGRERR_FAILURE) {
      datum.set_semi_major_axis(semi_major);
      datum.set_semi_minor_axis(semi_minor);
    }
    datum.meridian_offset() = gdal_spatial_ref.GetPrimeMeridian();
    // Set the proj4 string for datum.
    std::stringstream datum_proj4_ss;
    for(unsigned i=0; i < datum_strings.size(); i++)
      datum_proj4_ss << datum_strings[i] << ' ';
    // Add the current proj4 string in the case that our ellipse/datum
    // values are empty.
    if(boost::trim_copy(datum_proj4_ss.str()) == "")
      datum_proj4_ss << datum.proj4_str();
    datum.proj4_str() = boost::trim_copy(datum_proj4_ss.str());
    set_datum(datum);
  }
#endif // VW_HAVE_PKG_GDAL

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

    XY projected;
    LP unprojected;

    // Proj.4 expects the (lon,lat) pair to be in radians
    unprojected.u = lon_lat[0] * DEG_TO_RAD;
    unprojected.v = lon_lat[1] * DEG_TO_RAD;

    // Clamp the latitude range to [-HALFPI,HALFPI] ([-90, 90]) as occasionally
    // we get edge pixels that extend slightly beyond that range (probably due
    // to pixel as area vs point) and cause Proj.4 to fail. We use HALFPI
    // rather than other incantations for pi/2 because that's what proj.4 uses.
    if(unprojected.v > HALFPI)       unprojected.v = HALFPI;
    else if(unprojected.v < -HALFPI) unprojected.v = -HALFPI;

    projected = pj_fwd(unprojected, m_proj_context->proj_ptr());
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
    num_strings = boost::numeric_cast<int>(arg_strings.size());
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
}} // vw::cartography
