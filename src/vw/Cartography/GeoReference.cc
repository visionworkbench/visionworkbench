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
#include <boost/foreach.hpp>

// Proj.4
#include <proj_api.h>

// Macro for checking Proj.4 output, something we do a lot of.
#define CHECK_PROJ_ERROR(ctx_input) if(ctx_input.error_no()) vw_throw(ProjectionErr() << "Proj.4 error: " << pj_strerrno(ctx_input.error_no()))

namespace vw {
namespace cartography {

  bool read_georeference( GeoReference& georef,
                          ImageResource const& resource ) {

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

  void write_georeference( ImageResource& resource,
                           GeoReference const& georef ) {
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL* gdal =
      dynamic_cast<DiskImageResourceGDAL*>( &resource );
    if ( gdal ) return write_gdal_georeference( *gdal, georef );
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother
    // checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  bool read_header_string( ImageResource const& resource, std::string const& str_name,
                           std::string & str_val ) {

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL const* gdal =
      dynamic_cast<DiskImageResourceGDAL const*>( &resource );
    if ( gdal ) return read_gdal_string( *gdal, str_name, str_val );
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother
    // checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  void write_header_string( ImageResource& resource, std::string const& str_name,
                            std::string const& str_val ) {

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL* gdal =
      dynamic_cast<DiskImageResourceGDAL*>( &resource );
    if ( gdal ) write_gdal_string( *gdal, str_name, str_val );
    return;
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother
    // checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  std::string GeoReference::proj4_str() const {
    return m_proj_projection_str;
  }

  std::string GeoReference::overall_proj4_str() const {
    std::string proj4_str =
      m_proj_projection_str + " " + m_datum.proj4_str() + " +no_defs";
    return proj4_str;
  }

  void GeoReference::init_proj() {
    m_proj_context = ProjContext( overall_proj4_str() );
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
  GeoReference::GeoReference(Datum const& datum,
                             Matrix<double,3,3> const& transform) : GeoReferenceBase(datum) {
    set_transform(transform);
    set_geographic();
    init_proj();
  }

  GeoReference::GeoReference(Datum const& datum,
                             Matrix<double,3,3> const& transform,
                             PixelInterpretation pixel_interpretation) :
    GeoReferenceBase(datum, pixel_interpretation) {
    set_transform(transform);
    set_geographic();
    init_proj();
  }

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
  GeoReference::GeoReference(GeoReferenceDesc const& desc) {
    VW_ASSERT(desc.transform_size() == 9,
              IOErr() << "GeoReference::GeoReference: Unexpected number of elements in transform");

    m_datum = Datum(desc.datum());
    m_pixel_interpretation =
      static_cast<GeoReferenceBase::PixelInterpretation>(desc.pixel_interpretation());
    set_transform(Matrix3x3(desc.transform().data()));
    //m_is_projected = desc.is_projected();
    //m_proj_projection_str = desc.proj_projection_str();

    //init_proj();

    set_proj4_projection_str(desc.proj_projection_str());
  }

  GeoReferenceDesc GeoReference::build_desc() {
    GeoReferenceDesc desc;

    *(desc.mutable_datum()) = m_datum.build_desc();
    desc.set_pixel_interpretation(static_cast<GeoReferenceDesc::PixelInterpretation>(m_pixel_interpretation));
    std::copy(m_transform.begin(), m_transform.end(), RepeatedFieldBackInserter(desc.mutable_transform()));
    desc.set_is_projected(m_is_projected);
    desc.set_proj_projection_str(m_proj_projection_str);

    return desc;
  }
#endif

  void GeoReference::set_transform(Matrix3x3 transform) {
    m_transform = transform;
    m_shifted_transform = m_transform;
    m_shifted_transform(0,2) += 0.5*m_transform(0,0);
    m_shifted_transform(1,2) += 0.5*m_transform(1,1);
    m_inv_transform = vw::math::inverse(m_transform);
    m_inv_shifted_transform = vw::math::inverse(m_shifted_transform);

    update_lon_wrap();
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
    set_proj4_projection_str("+proj=longlat");
  }

  void GeoReference::set_equirectangular(double center_latitude, double center_longitude, double latitude_of_true_scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=eqc +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +lat_ts=" << latitude_of_true_scale
         << " +x_0=" << false_easting << " +y_0=" << false_northing
         << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_sinusoidal(double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sinu +lon_0=" << center_longitude << " +x_0="
         << false_easting << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_mercator(double center_latitude, double center_longitude, double latitude_of_true_scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=merc +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +lat_ts=" << latitude_of_true_scale
         << " +x_0=" << false_easting << " +y_0=" << false_northing
         << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=tmerc +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +k=" << scale << " +x_0=" << false_easting
         << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_orthographic(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=ortho +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +x_0=" << false_easting << " +y_0="
         << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_oblique_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sterea +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +k=" << scale << " +x_0=" << false_easting
         << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=stere +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +k=" << scale << " +x_0=" << false_easting
         << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=laea +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +x_0=" << false_easting << " +y_0="
         << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_lambert_conformal(double std_parallel_1, double std_parallel_2, double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=lcc +lat_1=" << std_parallel_1 << " +lat_2="
         << std_parallel_2 << " +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +x_0=" << false_easting << " +y_0="
         << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_UTM(int zone, int north) {
    std::ostringstream strm;
    strm << "+proj=utm +zone=" << zone;
    if (!north) strm << " +south";
    strm << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_proj4_projection_str(std::string const& s) {

    m_proj_projection_str = s; // Store the string in this class (it is also stored in m_proj_context)

    // Extract some information from the string
    if (s.find("+proj=longlat") == 0)
      m_is_projected = false;
    else
      m_is_projected = true;

    // Force all "eqc" projections to turn off -180 to 180 longitude wrapping in proj4.
    //// - This is only temporary so we can make another modification below
    if ( (s.find("+proj=eqc") == 0) && (s.find("+over") == std::string::npos) )
      m_proj_projection_str.append(" +over");

    init_proj(); // Initialize m_proj_context

    update_lon_wrap();
  }


  void GeoReference::update_lon_wrap() {
  
      if (m_proj_projection_str.find("+proj=eqc") != 0) {
        m_using_lon_wrap = false; // Other projections currently not using this correction.
        return; // Nothing else to do here!
     }
     
     // Start of special handling code for eqc case.
     
     // Since we don't have any image information we have to assume this
     //  georef needs to be valid for a full 360 degrees.
     const double GEOREF_VALID_WIDTH = 360.0;
     const double halfWidth = GEOREF_VALID_WIDTH / 2.0;

     // Proj4 won't work with angles outside of these.  Pulled from pj_fwd.cc
     const double PROJ4_MAX_LON =  10 * RAD_TO_DEG;
     const double PROJ4_MIN_LON = -10 * RAD_TO_DEG;

     // This method won't work if there is a rotation so check for that.
     if (fabs(m_transform(0,1)) > 0.01)
       vw_throw(NoImplErr() << "EQC projections with rotation are not supported.");

     // Figure out where the 0,0 pixel transforms to in lon/lat.
     // - Remember that we turned off wrapping before this call.
     Vector2 lonlatBound = pixel_to_lonlat(Vector2(0,0)); // No wiggle room here since we are aligning to 90 degrees anyways
     double minLon, maxLon;                                          

     // In order to make the normalization range look better, shift it
     //  to start at the nearest multiple of 90 degrees.
     // - Doing this will cause us to fail on 360 degree images that do not start on a multiple of 90!
     const double ALIGN_MULTIPLE = 90.0;
     if (m_transform(0,0) > 0) { // This is the minimum longitude
       minLon = GEOREF_VALID_WIDTH * 4; // Start outside of the valid Proj4 bounds
       while (minLon > lonlatBound[0]) // Get the next multiple of 90 below the bound we found
         minLon -= ALIGN_MULTIPLE;
       maxLon = minLon + GEOREF_VALID_WIDTH;
     }
     else { // This is the maximum longitude, offset to get the minimum.
       maxLon = -GEOREF_VALID_WIDTH * 4; // Start outside of the valid Proj4 bounds
       while (maxLon < lonlatBound[0]) // Get the next multiple of 90 below the bound we found
         maxLon += ALIGN_MULTIPLE;
       minLon = maxLon - GEOREF_VALID_WIDTH;
     }  
     
     // Now that we now the range that the the georef "naturally"
     //  projects from, get the center point.

     // Need to adjust to make sure we stay inside proj4 bounds
     // - A better solution is needed to function outside those bounds, but that would add a lot of complexity.
     if (minLon < PROJ4_MIN_LON) { // Shift upwards to -10 radians
       minLon = PROJ4_MIN_LON;
     }
     if (maxLon > PROJ4_MAX_LON) { // Shift downwards to 10 radians
       minLon -= (maxLon - PROJ4_MAX_LON);
     }

     m_center_lon_wrap = minLon + halfWidth;
     m_using_lon_wrap  = true;
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

    std::vector<std::string> input_strings, output_strings, datum_strings;
    std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
    boost::split( input_strings, trimmed_proj4_str, boost::is_any_of(" ") );
    for (size_t i = 0; i < input_strings.size(); ++i) {
      const std::string& key = input_strings[i];

      // Pick out the parts of the projection string that pertain to
      // map projections.  We essentially want to eliminate all of
      // the strings that have to do with the datum, since those are
      // handled by interacting directly with the
      // OGRSpatialReference below. This is sort of messy, but it's
      // the easiest way to do this, as far as I can tell.
      if (key == "+k=0") {
        vw_out(WarningMessage) << "Input contained an illegal scale_factor of zero. Ignored." << std::endl;
      } else if ((key.find("+proj=") == 0) ||
          (key.find("+x_0=") == 0) ||
          (key.find("+y_0=") == 0) ||
          (key.find("+lon") == 0) ||
          (key.find("+lat") == 0) ||
          (key.find("+k=") == 0) ||
          (key.find("+lat_ts=") == 0) ||
          (key.find("+ns") == 0) ||
          (key.find("+no_cut") == 0) ||
          (key.find("+h=") == 0) ||
          (key.find("+W=") == 0) ||
          (key.find("+units=") == 0) ||
          (key.find("+zone=") == 0)) {
        output_strings.push_back(key);
      } else if ((key.find("+ellps=") == 0) ||
                 (key.find("+datum=") == 0)) {
        // We put these in the proj4_str for the Datum class.
        datum_strings.push_back(key);
      }
    }
    std::ostringstream strm;
    BOOST_FOREACH( std::string const& element, output_strings )
      strm << element << " ";

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
    BOOST_FOREACH( std::string const& element, datum_strings )
      datum_proj4_ss << element << " ";
    // Add the current proj4 string in the case that our ellipse/datum
    // values are empty.
    if ( boost::trim_copy(datum_proj4_ss.str()) == "" )
      datum_proj4_ss << datum.proj4_str();
    datum.proj4_str() = boost::trim_copy(datum_proj4_ss.str());
    set_datum(datum);
  }
#endif // VW_HAVE_PKG_GDAL

  /// For a given pixel coordinate, compute the position of that
  /// pixel in this georeferenced space.
  Vector2 GeoReference::pixel_to_point(Vector2 pix) const {
    Vector2 loc;
    Matrix3x3 M = this->vw_native_transform();
    double denom = pix[0] * M(2,0) + pix[1] * M(2,1) + M(2,2);
    loc[0] = (pix[0] * M(0,0) + pix[1] * M(0,1) + M(0,2)) / denom;
    loc[1] = (pix[0] * M(1,0) + pix[1] * M(1,1) + M(1,2)) / denom;
    return loc;
  }

  /// For a given location 'loc' in projected space, compute the
  /// corresponding pixel coordinates in the image.
  Vector2 GeoReference::point_to_pixel(Vector2 loc) const {
    Vector2 pix;
    Matrix3x3 M = this->vw_native_inverse_transform();
    double denom = loc[0] * M(2,0) + loc[1] * M(2,1) + M(2,2);
    pix[0] = (loc[0] * M(0,0) + loc[1] * M(0,1) + M(0,2)) / denom;
    pix[1] = (loc[0] * M(1,0) + loc[1] * M(1,1) + M(1,2)) / denom;
    return pix;
  }


  /// For a point in the projected space, compute the position of
  /// that point in unprojected (Geographic) coordinates (lat,lon).
  Vector2 GeoReference::point_to_lonlat(Vector2 loc) const {
    if ( ! m_is_projected ) return loc;

    projXY projected;
    projLP unprojected;

    projected.u = loc[0]; // Store in proj4 object
    projected.v = loc[1];

    // Call proj4 to do the conversion and check for errors.
    unprojected = pj_inv(projected, m_proj_context.proj_ptr());
    CHECK_PROJ_ERROR( m_proj_context );

    // Convert from radians to degrees.
    return Vector2(unprojected.u * RAD_TO_DEG, unprojected.v * RAD_TO_DEG);
  }


  /// Given a position in geographic coordinates (lat,lon), compute
  /// the location in the projected coordinate system.
  Vector2 GeoReference::lonlat_to_point(Vector2 lon_lat) const {
    if ( ! m_is_projected ) return lon_lat;

    // For eqc projections, transform intput lon into requested range.
    if (m_using_lon_wrap) {
      const double OFFSET = 360.0;
      const double RANGE  = 180.0;
      const double maxLon = m_center_lon_wrap + RANGE;
      const double minLon = m_center_lon_wrap - RANGE;

      while (lon_lat[0] > maxLon) // Too high
        lon_lat[0] -= OFFSET;
      while (lon_lat[0] < minLon) // Too low
        lon_lat[0] += OFFSET;
    }


    // This value is proj's internal limit
    static const double BOUND = 1.5707963267948966 - (1e-10) - std::numeric_limits<double>::epsilon();

    projXY projected;
    projLP unprojected;

    // Proj.4 expects the (lon,lat) pair to be in radians
    unprojected.u = lon_lat[0] * DEG_TO_RAD;
    unprojected.v = lon_lat[1] * DEG_TO_RAD;

    // Clamp the latitude range to [-HALFPI,HALFPI] ([-90, 90]) as occasionally
    // we get edge pixels that extend slightly beyond that range (probably due
    // to pixel as area vs point) and cause Proj.4 to fail. We use HALFPI
    // rather than other incantations for pi/2 because that's what proj.4 uses.
    if(unprojected.v > BOUND)        unprojected.v = BOUND;
    else if(unprojected.v < -BOUND) unprojected.v = -BOUND;

    // Call proj4 to do the conversion and check for errors.
    projected = pj_fwd(unprojected, m_proj_context.proj_ptr());
    CHECK_PROJ_ERROR( m_proj_context );

    return Vector2(projected.u, projected.v);
  }


  //*****************************************************************
  //************** Functions for class ProjContext ******************

  char** ProjContext::split_proj4_string(std::string const& proj4_str, int &num_strings) {
    std::vector<std::string> arg_strings;
    std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
    boost::split( arg_strings, trimmed_proj4_str, boost::is_any_of(" ") );

    char** strings = new char*[arg_strings.size()];
    for ( size_t i = 0; i < arg_strings.size(); ++i ) {
      strings[i] = new char[2048];
      strncpy(strings[i], arg_strings[i].c_str(), 2048);
    }
    num_strings = boost::numeric_cast<int>(arg_strings.size());
    return strings;
  }

#if PJ_VERSION < 480
  ProjContext::ProjContext(std::string const& proj4_str) : m_proj4_str(proj4_str) {

    // proj.4 is expecting the parameters to be split up into seperate
    // c-style strings.
    int num;
    char** proj_strings = split_proj4_string(m_proj4_str, num);
    m_proj_ptr.reset( pj_init(num, proj_strings),
                      pj_free );

    VW_ASSERT( !pj_errno, InputErr() << "Proj.4 failed to initialize on string: " << m_proj4_str << "\n\tError was: " << pj_strerrno(pj_errno) );

    for ( int i = 0; i < num; i++ ) delete [] proj_strings[i];
    delete [] proj_strings;
  }
  ProjContext::ProjContext( ProjContext const& other ) : m_proj_ptr(other.m_proj_ptr), m_proj4_str(other.m_proj4_str) {}
  int ProjContext::error_no() const {
    return pj_errno;
  }
#else // PJ_VERSION >= 480

  ProjContext::ProjContext(std::string const& proj4_str ) : m_proj4_str(proj4_str) {
    m_proj_ctx_ptr.reset(pj_ctx_alloc(),pj_ctx_free);
    int num;
    char** proj_strings = split_proj4_string(m_proj4_str, num);
    m_proj_ptr.reset(pj_init_ctx( m_proj_ctx_ptr.get(),
                                  num, proj_strings ),
                     pj_free);

    VW_ASSERT( !pj_ctx_get_errno(m_proj_ctx_ptr.get()),
               InputErr() << "Proj.4 failed to initialize on string: " << m_proj4_str << "\n\tError was: " << pj_strerrno(pj_ctx_get_errno(m_proj_ctx_ptr.get())) );

    for ( int i = 0; i < num; i++ ) delete [] proj_strings[i];
    delete [] proj_strings;
  }

  ProjContext::ProjContext( ProjContext const& other ) : m_proj4_str(other.m_proj4_str) {
    m_proj_ctx_ptr.reset(pj_ctx_alloc(),pj_ctx_free);
    if ( m_proj4_str.empty() )
      return; // They've made a copy of an uninitialized
              // projcontext. Not an error .. since they can
              // initialize later.

    int num;
    char** proj_strings = split_proj4_string(m_proj4_str, num);
    m_proj_ptr.reset(pj_init_ctx( m_proj_ctx_ptr.get(),
                                  num, proj_strings ),
                     pj_free);

    VW_ASSERT( !pj_ctx_get_errno(m_proj_ctx_ptr.get()),
               InputErr() << "Proj.4 failed to initialize on string: " << m_proj4_str << "\n\tError was: " << pj_strerrno(pj_ctx_get_errno(m_proj_ctx_ptr.get())) );

    for ( int i = 0; i < num; i++ ) delete [] proj_strings[i];
    delete [] proj_strings;
  }

  int ProjContext::error_no() const {
    return pj_ctx_get_errno(m_proj_ctx_ptr.get());
  }
#endif
//************** End functions for class ProjContext ******************
//*********************************************************************


  // Simple GeoReference modification tools


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

  GeoReference resample( GeoReference const& input,
                         double scale_x, double scale_y ) {
    GeoReference output = input;
    Matrix3x3 T = output.transform();
    T(0,0) /= scale_x;
    T(1,1) /= scale_y;
    if ( input.pixel_interpretation() == GeoReference::PixelAsArea ) {
      Vector2 top_left_ll =
        input.pixel_to_point( -Vector2(0.5 / scale_x, 0.5 / scale_y) );
      T(0,2) = top_left_ll[0];
      T(1,2) = top_left_ll[1];
    }
    output.set_transform(T);
    return output;
  }

  GeoReference resample( GeoReference const& input,
                         double scale ) {
    return resample(input, scale, scale );
  }

}} // vw::cartography

#undef CHECK_PROJ_ERROR
