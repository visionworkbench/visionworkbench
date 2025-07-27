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
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/FileTypes.h>
#include <vw/Core/StringUtils.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/BresenhamLine.h>

#include <ogr_spatialref.h>
#include <cpl_string.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

// TODO(oalexan1): Wipe all mention of proj4_str.
// Everything must be done by setting the OGRSpatialReference object.
// TODO(oalexan1): Look at all leftover todo statements.
// TODO(oalexan1): Look in GeoTransform for things to fix.

// Proj
// TODO(oalexan1): Wipe proj code, now the gdal wrapper is used

// TODO(oalexan1): Wipe from everywhere the invocation of VW_HAVE_PKG_GDAL
// VW will always be built with GDAL support.

namespace vw {
namespace cartography {

// Update a georeference based on an srs string and/or datum.
// This function is more likely to remember the datum name than set_wkt().
void set_srs_string(std::string srs_string, bool have_user_datum,
                    vw::cartography::Datum const& user_datum,
                    vw::cartography::GeoReference & georef) {

  // When an EPSG code is provided, store the name so that
  //  it shows up when the GeoReference object is written
  //  out to disk.
  if (srs_string.find("EPSG") != std::string::npos)
    georef.set_projcs_name(srs_string);

  // Set srs_string into given georef. Note that this may leave the
  // georef's affine transform inconsistent.

  // TODO: The line below needs more thought
  if (srs_string == "")
    srs_string = "+proj=longlat";

  // TODO(oalexan1): It is not clear this will be enough to remember the datum
  // name. May need to call set_datum() at the end, if have_user_datum is true,
  // and if the georef does not have the datum name but the datum has it.
  if (have_user_datum)
    srs_string += " " + user_datum.proj4_str();
  
  OGRSpatialReference gdal_spatial_ref;
  if (gdal_spatial_ref.SetFromUserInput(srs_string.c_str()) != OGRERR_NONE)
    vw::vw_throw(vw::ArgumentErr() << "Failed to parse: \"" << srs_string << "\".");
  char *wkt_str_tmp = NULL;
  gdal_spatial_ref.exportToWkt(&wkt_str_tmp);
  srs_string = wkt_str_tmp;
  CPLFree(wkt_str_tmp);
  georef.set_wkt(srs_string);
}

bool read_georeference(GeoReference& georef,
                        ImageResource const& resource) {
  DiskImageResourceGDAL const* gdal =
    dynamic_cast<DiskImageResourceGDAL const*>(&resource);
  if (gdal) 
    return read_gdal_georeference(georef, *gdal);

  DiskImageResourcePDS const* pds =
    dynamic_cast<DiskImageResourcePDS const*>(&resource);
  if(pds) 
    return read_pds_georeference(georef, *pds);
  return false;
}

/// A convenience function to read georeferencing information from an image file.
bool read_georeference(GeoReference& georef, const std::string &filename) {
  
  // No image with a SPOT5 suffix can ever have georeference.
  if (vw::has_spot5_extension(filename)) 
    return false;

  boost::shared_ptr<DiskImageResource> r(DiskImageResourcePtr(filename));
  bool result = read_georeference(georef, *r);
  
  return result;
}

void write_georeference(ImageResource& resource,
                          GeoReference const& georef) {
  DiskImageResourceGDAL* gdal = dynamic_cast<DiskImageResourceGDAL*>(&resource);
  if (gdal) 
    return write_gdal_georeference(*gdal, georef);

  // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
  vw_throw(NoImplErr()
            << "This image resource does not support writing georeferencing information.");
}

bool read_header_string(ImageResource const& resource, std::string const& str_name,
                          std::string & str_val) {

  DiskImageResourceGDAL const* gdal = dynamic_cast<DiskImageResourceGDAL const*>(&resource);
  if (gdal) 
    return read_gdal_string(*gdal, str_name, str_val);

  // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
  vw_throw(NoImplErr()
            << "This image resource does not support writing georeferencing information.");
}

/// A function to read the value of a variable with given name from a geotiff
/// file header. Return an empty string on failure.
std::string read_header_string(std::string filename, std::string const& str_name) {
  boost::shared_ptr<DiskImageResource> r(DiskImageResourcePtr(filename));
  std::string str_val;
  read_header_string(*r, str_name, str_val);
  return str_val;
}

bool read_header_strings(ImageResource const& resource, 
                          std::map<std::string, std::string> & value_pairs) {

  DiskImageResourceGDAL const* gdal = dynamic_cast<DiskImageResourceGDAL const*>(&resource);
  if (gdal) 
    return read_gdal_strings(*gdal, value_pairs);
  // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
  vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
}

void write_header_string(ImageResource& resource, std::string const& str_name,
                          std::string const& str_val) {
  DiskImageResourceGDAL* gdal =
    dynamic_cast<DiskImageResourceGDAL*>(&resource);
  if (gdal) write_gdal_string(*gdal, str_name, str_val);
  return;
  // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
  vw_throw(NoImplErr()
            << "This image resource does not support writing georeferencing information.");
}

// Remove duplicates after concatenation of proj4 strings. Keep the first
// instance, as the second comes from the datum string, which may not be right.

size_t remove_proj4_duplicates(std::string const& str_in, std::string &str_out) {
  std::vector<std::string> arg_strings;
  std::string trimmed_proj4_str = boost::trim_copy(str_in);
  boost::split(arg_strings, trimmed_proj4_str, boost::is_any_of(" "));
  str_out = "";

  std::set<int> duplicates;
  for (size_t i = 0; i < arg_strings.size(); i++) {

    bool duplicate = false;
    for (size_t k = i + 1; k < arg_strings.size(); k++) {
      if (arg_strings[i] == arg_strings[k]) {
        duplicates.insert(k);
        break;
      }

      // Do not let the projection show up twice. The second one is the default one,
      // and may not be right.
      if (arg_strings[i].rfind("+proj=", 0) == 0  && 
          arg_strings[k].rfind("+proj=", 0) == 0) {
        duplicates.insert(k);
        break;
      }
    
    } // End k loop
  }

  size_t num_kept = 0;
  for (size_t i = 0; i < arg_strings.size(); i++) {
  
    if (duplicates.find(i) != duplicates.end())
      continue; // duplicate
    
    if (num_kept > 0)
      str_out += " ";
  
    str_out += arg_strings[i];
    num_kept++;
  } // End i loop

  return num_kept;
}
  
// If there are no names for the datum and ellipsoid, copy
// from this datum, if the params agree. Note that WGS84
// and GRS 80 (NAD83) semi-major axes are the same, but the
// semi-minor axes differ by a bit more than 1e-4.
// Return true if the data copy took place.
bool copyFromDatum(OGRSpatialReference & gdal_spatial_ref, 
                   vw::cartography::Datum const& datum) {

  // The datum name can be null, and this can result in a crash.
  // Cannot continue in this case.
  const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
  if (datum_name == NULL)
    return false;
  
  std::string lc_datum_name = boost::to_lower_copy(std::string(datum_name));
  if (lc_datum_name.find("unknown") == std::string::npos &&
      lc_datum_name.find("user specified datum") == std::string::npos)
    return false; // The datum name is already good
    
  const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
  if (spheroid_name == NULL)
    return false;
  const char* meridian_name_ptr = gdal_spatial_ref.GetAttrValue("PRIMEM");
  if (meridian_name_ptr == NULL)
    return false;
  std::string meridian_name = meridian_name_ptr;
  double meridian_offset = gdal_spatial_ref.GetPrimeMeridian();
  double inv_flattening = gdal_spatial_ref.GetInvFlattening();
  OGRErr e1, e2;
  double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
  double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
  if (e1 == OGRERR_FAILURE || e2 == OGRERR_FAILURE)
    vw::vw_throw(vw::ArgumentErr() << "Failed to red axes from OGRSpatialRef.\n");

  if (datum.semi_major_axis() == semi_major &&
      std::abs(datum.semi_minor_axis() - semi_minor) < 1e-6 &&
      boost::to_lower_copy(datum.meridian_name()) == boost::to_lower_copy(meridian_name) &&
      datum.meridian_offset() == meridian_offset) {
    gdal_spatial_ref.SetGeogCS("Geographic Coordinate System", 
                                datum.name().c_str(),
                                datum.spheroid_name().c_str(),
                                semi_major,
                                inv_flattening,
                                datum.meridian_name().c_str(),
                                meridian_offset);
     return true;
  }
  
   return false;
}

// Try to fix datum name in the spatial ref by comparing to several
// known datums.
void fixDatum(OGRSpatialReference & gdal_spatial_ref) {
  
  // The datum name can be null, and this can result in a crash.
  const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
  if (datum_name == NULL) 
    return;

  std::string lc_datum_name = boost::to_lower_copy(std::string(datum_name));
  if (lc_datum_name.find("unknown") == std::string::npos &&
      lc_datum_name.find("user specified datum") == std::string::npos)
    return; // The datum name is already good

  OGRErr e;
  double semi_major = gdal_spatial_ref.GetSemiMajor(&e);
  if (e == OGRERR_FAILURE)
    vw::vw_throw(vw::ArgumentErr() << "Failed to red axes from OGRSpatialRef.\n");
    
  if (semi_major == 6378137) {
    // Try WGS84
    vw::cartography::Datum datum;
    datum.set_well_known_datum("WGS84");
    copyFromDatum(gdal_spatial_ref, datum);
  } else if (semi_major == 1737400) {
    // Try the moon
    vw::cartography::Datum datum;
    datum.set_well_known_datum("D_MOON");
    copyFromDatum(gdal_spatial_ref, datum);
  } else if (semi_major == 3396190) {
    // Try Mars
    vw::cartography::Datum datum;
    datum.set_well_known_datum("D_MARS");
    copyFromDatum(gdal_spatial_ref, datum);
  } 

  return;  
}

// Low-level function
std::string ogr_wkt(OGRSpatialReference const & ogr) {
  char *wkt;
  ogr.exportToWkt(&wkt);
  std::string result(wkt);
  CPLFree(wkt);
  return result;
}

std::string GeoReference::proj4_str() const {
  return m_proj_projection_str;
}

std::string GeoReference::overall_proj4_str() const {
  // Make sure these elements exist but prevent duplicate entries
  std::string proj4_str = boost::trim_copy(m_proj_projection_str) + " "
                          + boost::trim_copy(m_datum.proj4_str()) + " +no_defs";
  std::string proj4_str_no_dups;
  remove_proj4_duplicates(proj4_str, proj4_str_no_dups);

  return proj4_str_no_dups;
}

// This will recreate the GeoReference object
void GeoReference::init_proj() {
  // Update the projection context object with the current proj4 string, 
  //  then make sure the lon center is still correct.
    
  // This will append the datum info (but not the name)
  std::string srs_string = overall_proj4_str(); 
  
  // Form the georeference from the proj4 string  
  if (m_gdal_spatial_ref.SetFromUserInput(srs_string.c_str()) != OGRERR_NONE)
    vw::vw_throw(vw::ArgumentErr() << "Failed to parse: " << srs_string << "\n");
  set_wkt(ogr_wkt(m_gdal_spatial_ref));
}

// The empty constructor. This initialises m_gdal_spatial_ref to an empty
// object, and m_gdal_spatial_ref.exportToWkt() will return an empty string.
GeoReference::GeoReference(): m_pixel_interpretation(PixelAsArea), m_projcs_name("") {
  
  set_transform(vw::math::identity_matrix<3>());
  set_geographic(); // will call init_proj()
}

GeoReference::GeoReference(Datum const& datum):
      m_pixel_interpretation(PixelAsArea), m_datum(datum), m_projcs_name("") {
  set_transform(vw::math::identity_matrix<3>());
  set_geographic(); // will call init_proj()
}

GeoReference::GeoReference(Datum const& datum, PixelInterpretation pixel_interpretation):
 m_pixel_interpretation (pixel_interpretation), m_datum(datum), m_projcs_name("") {
  set_transform(vw::math::identity_matrix<3>());
  set_geographic(); // will call init_proj()
}

GeoReference::GeoReference(Datum const& datum,
                            Matrix<double,3,3> const& transform):
                  m_pixel_interpretation(PixelAsArea), m_datum(datum), m_projcs_name("") {
  set_transform(transform);
  set_geographic(); // will call init_proj()
}

GeoReference::GeoReference(Datum const& datum,
                            Matrix<double,3,3> const& transform,
                            PixelInterpretation pixel_interpretation):
  m_pixel_interpretation(pixel_interpretation), m_datum(datum), m_projcs_name("") {
  set_transform(transform);
  set_geographic(); // will call init_proj()
}

// Destructor
GeoReference::~GeoReference() {
}

void GeoReference::set_transform(Matrix3x3 transform) {
  m_transform = transform;
  m_shifted_transform = m_transform;
  m_shifted_transform(0,2) += 0.5*m_transform(0,0);
  m_shifted_transform(1,2) += 0.5*m_transform(1,1);
  m_inv_transform         = vw::math::inverse(m_transform);
  m_inv_shifted_transform = vw::math::inverse(m_shifted_transform);

  // If the georef is already fully initialized, update lon-lat box.
  // Otherwise this is called from a constructor, and a subsequent 
  // step will result in this being set.
  if (m_proj_context.is_initialized())
    ll_box_from_pix_box(vw::BBox2(0, 0, 2, 2)); // must be at least 2 pixels in size
}

// Set the datum. Keep the projection.
void GeoReference::set_datum(Datum const& datum) {
  
  m_datum = datum;
  std::string projcs_name = m_projcs_name; 
  if (projcs_name.empty())
    projcs_name = "Geographic Coordinate System";

  m_gdal_spatial_ref.SetGeogCS(projcs_name.c_str(), 
                                datum.name().c_str(),
                                datum.spheroid_name().c_str(),
                                datum.semi_major_axis(),
                                datum.inverse_flattening(),
                                datum.meridian_name().c_str(),
                                datum.meridian_offset());
  
  // Recreate the georeference
  set_wkt(ogr_wkt(m_gdal_spatial_ref));
}

// Adjust the affine transform to the VW convention ([0,0] is at
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
  this->set_datum(m_datum); // set the datum, this will rebuild the georef
}

void GeoReference::set_geographic() {
  set_proj4_projection_str("+proj=longlat");
}

void GeoReference::set_equirectangular(double center_latitude, double center_longitude,
                                        double latitude_of_true_scale, double false_easting,
                                        double false_northing) {
  std::ostringstream strm;
  strm << "+proj=eqc +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +lat_ts=" << latitude_of_true_scale
        << " +x_0=" << false_easting << " +y_0=" << false_northing
        << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_sinusoidal(double center_longitude, double false_easting,
                                  double false_northing) {
  std::ostringstream strm;
  strm << "+proj=sinu +lon_0=" << center_longitude << " +x_0="
        << false_easting << " +y_0=" << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_mercator(double center_latitude, double center_longitude,
                                double latitude_of_true_scale, double false_easting,
                                double false_northing) {
  std::ostringstream strm;
  strm << "+proj=merc +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +lat_ts=" << latitude_of_true_scale
        << " +x_0=" << false_easting << " +y_0=" << false_northing
        << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_transverse_mercator(double center_latitude, double center_longitude,
                                            double scale, double false_easting,
                                            double false_northing) {
  std::ostringstream strm;
  strm << "+proj=tmerc +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +k=" << scale << " +x_0=" << false_easting
        << " +y_0=" << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_orthographic(double center_latitude, double center_longitude,
                                    double false_easting, double false_northing) {
  std::ostringstream strm;
  strm << "+proj=ortho +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +x_0=" << false_easting << " +y_0="
        << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_stereographic(double center_latitude, double center_longitude,
                                      double scale, double false_easting, double false_northing) {
  std::ostringstream strm;
  strm << "+proj=stere +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +k=" << scale << " +x_0=" << false_easting
        << " +y_0=" << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_oblique_stereographic(double center_latitude,
                                              double center_longitude, double scale,
                                              double false_easting, double false_northing) {
  std::ostringstream strm;
  strm << "+proj=sterea +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +k=" << scale << " +x_0=" << false_easting
        << " +y_0=" << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_gnomonic(double center_latitude, double center_longitude, 
                                double scale,
                                double false_easting, double false_northing) {
  
  std::ostringstream strm;
  strm << "+proj=gnom +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +k=" << scale << " +x_0=" << false_easting
        << " +y_0=" << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_lambert_azimuthal(double center_latitude, double center_longitude,
                                          double false_easting, double false_northing) {
  std::ostringstream strm;
  strm << "+proj=laea +lon_0=" << center_longitude << " +lat_0="
        << center_latitude << " +x_0=" << false_easting << " +y_0="
        << false_northing << " +units=m";
  set_proj4_projection_str(strm.str());
}

void GeoReference::set_lambert_conformal(double std_parallel_1, double std_parallel_2,
                                          double center_latitude, double center_longitude,
                                          double false_easting, double false_northing) {
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

// This function does not set the datum radius based on the projection, which is wrong.
// Use instead set_srs_string().
// TODO(oalexan1): This must be made private as it does not fully create
// the object.
void GeoReference::set_proj4_projection_str(std::string const& s) {
  
  // Store the string in this class
  m_proj_projection_str = boost::trim_copy(s);

  // Extract some information from the string
  // TODO(oalexan1): This is duplicate text that must be removed
  if (m_proj_projection_str.find("+proj=longlat") == 0)
    m_is_projected = false;
  else
    m_is_projected = true;

  // TODO(oalexan1): Wipe this code
  // Disable -180 to 180 longitude wrapping in proj4.
  // - With wrapping off, Proj4 can work significantly outside those ranges (though there is a limit)
  // - We will make sure that the input longitudes are in a safe range.
  if  ((m_proj_projection_str.find("+over") == std::string::npos) &&
        (m_proj_projection_str.find("+proj=utm") == std::string::npos))
    m_proj_projection_str.append(" +over");

  init_proj(); // Initialize the projection
}

// TODO(oalexan1): Wipe this
bool GeoReference::extract_proj4_value(std::string const& proj4_string, 
                                       std::string const& key,
                                       std::string &s) {
  // Try to find the key
  size_t key_pos = proj4_string.find(key);
  if (key_pos == std::string::npos)
    return false;
  size_t key_end = key_pos + key.size();

  // Figure out the bounds of the number
  size_t eq_pos    = proj4_string.find("=", key_pos);
  size_t space_pos = proj4_string.find(" ", eq_pos);
  if ((eq_pos == std::string::npos) ||
        (eq_pos - key_end > 2)) // Make sure we got the right "="
    return false;
  if (space_pos == std::string::npos)
    space_pos = proj4_string.size();
  size_t start  = eq_pos + 1;
  size_t length = space_pos - start;
  s = proj4_string.substr(eq_pos+1, length);
  return true;
}

// TODO(oalexan1): Wipe this
bool GeoReference::extract_proj4_value(std::string const& proj4_string, 
                                       std::string const& key,
                                       double &value) {
  std::string s;
  if (!extract_proj4_value(proj4_string, key, s))
    return false;
  value = atof(s.c_str());
  return true;
}

// Strip the "+over" text from our stored proj4 info, but don't ll_box_from_pix_box().
// - Used to strip an extra tag out of [-180,180] range images where it is not needed.
void GeoReference::clear_proj4_over() {
  
  // This requires an initialized georeference  
  if (!m_proj_context.is_initialized())
    vw_throw(NoImplErr() << "GeoReference::clear_proj4_over() requires an initialized georeference.");
    
  return;
  // TODO(oalexan1): Remove all dependencies on proj4_str.
  // Clear out m_proj_projection_str, then recreate the ProjContext object.
  if (string_replace(m_proj_projection_str, "+over", "")) {
    // If we had to make any changes, strip out any double spaces and 
    //  trailing spaces and then update our ProjContext object.
    string_replace(m_proj_projection_str, "  ", " ");
    m_proj_projection_str = boost::trim_copy(m_proj_projection_str);
    init_proj();
  }
}

// Add the "+over" text to our stored proj4 info, but don't ll_box_from_pix_box().
void GeoReference::set_proj4_over() {
  return;
  // TODO(oalexan1): Wipe this
  // Clear out m_proj_projection_str, then recreate the ProjContext object.
  if (m_proj_projection_str.find("+over") == std::string::npos) {
    m_proj_projection_str.append(" +over");
  
    // If we had to make any changes, strip out any double spaces and 
    //  trailing spaces and then update our ProjContext object.
    string_replace(m_proj_projection_str, "  ", " ");
    m_proj_projection_str = boost::trim_copy(m_proj_projection_str);
    //m_proj_context = ProjContext(overall_proj4_str(), geo_wkt);
    init_proj();
  }
}

double GeoReference::test_pixel_reprojection_error(Vector2 const& pixel) {
  Vector2 out_pixel = lonlat_to_pixel(pixel_to_lonlat(pixel));
  Vector2 diff = out_pixel - pixel;
  double error = sqrt(diff.x()*diff.x() + diff.y()*diff.y());
  return error;
}

// Every function that modifies the georef must call this function
void GeoReference::set_wkt(std::string const& wkt) {
  
  if (wkt.empty())
    return;
  
  // A copy of the existing datum, before the reset.
  vw::cartography::Datum prior_datum = m_datum;

  // Set the spatial reference object from the WKT string  
  m_gdal_spatial_ref.importFromWkt(wkt.c_str());
  
  // If the datum name is not known, try to set it based on a datum with the same
  // parameters.
  fixDatum(m_gdal_spatial_ref);

  // If there is a PROJCS name, record it.
  m_projcs_name = "Geographic Coordinate System";
  const char * projcs = m_gdal_spatial_ref.GetAttrValue("PROJCS");
  if (projcs != NULL)
    m_projcs_name = std::string(projcs); // Careful here, to avoid a segfault

  // The returned coordinates will be in longitude, latitude order
  // https://gdal.org/tutorials/osr_api_tut.html#coordinate-transformation
  m_gdal_spatial_ref.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  
  // Create a new ProjContext object. This one has pointers to transforms, so a
  // new object is created to force a clean start. Same will happen on copy.
  // Create the underlying datum. Use a local scope to force a deallocation.
  m_proj_context = ProjContext(); 
  m_proj_context.m_proj_crs = m_gdal_spatial_ref; // need a copy of this
  {
    boost::shared_ptr<OGRSpatialReference> geoCS(m_gdal_spatial_ref.CloneGeogCS()); 
    geoCS->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    m_datum.set_datum_from_spatial_ref(*geoCS);
    m_proj_context.m_lonlat_crs = *geoCS; // need a copy of this
  }
  
  // This is a fix for GDAL forgetting datum names for small bodies when setting
  // a new projection.
  bool copied = copyFromDatum(m_gdal_spatial_ref, prior_datum);
  if (copied)
    m_datum.set_datum_from_spatial_ref(m_gdal_spatial_ref);
  
  // Create the lonlat vs proj coordinate transformations
  m_proj_context.init_transforms();  

  // Find the lon-lat bbox. This will be updated later if the transform
  // is set or if the georef is read from disk.
  // TODO(oalexan1): Record the image width and height on loading and use here.
  ll_box_from_pix_box(vw::BBox2(0, 0, 2, 2)); // must be at least 2 pixels in size
  
  // TODO(oalexan1): May need to wipe all below. Now that everything is based off
  // the spatial ref, there may be no need for this. 
  
  // Read projection information out of the spatial ref. TODO(oalexan1): Wipe
  // this. The ultimate source of truth is the spatial ref.
  char* proj4_str_tmp;
  m_gdal_spatial_ref.exportToProj4(&proj4_str_tmp);
  std::string proj4_str = proj4_str_tmp;
  CPLFree(proj4_str_tmp);

  std::vector<std::string> input_strings, output_strings, datum_strings;
  std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
  boost::split(input_strings, trimmed_proj4_str, boost::is_any_of(" "));
  BOOST_FOREACH(const std::string& key, input_strings) {
    // Pick out the parts of the projection string that pertain to
    // map projections.  We essentially want to eliminate all of
    // the strings that have to do with the datum, since those are
    // handled by interacting directly with the
    // OGRSpatialReference below. This is sort of messy, but it's
    // the easiest way to do this, as far as I can tell.
    if (key == "+k=0") {
      vw_out(WarningMessage)
          << "Input contained an illegal scale_factor of zero. Ignored."
          << std::endl;
    } else if (boost::starts_with(key, "+proj=") ||
                boost::starts_with(key, "+x_0=") ||
                boost::starts_with(key, "+y_0=") ||
                boost::starts_with(key, "+lon") ||
                boost::starts_with(key, "+lat") ||
                boost::starts_with(key, "+k=") ||
                boost::starts_with(key, "+k_0=") ||
                boost::starts_with(key, "+lat_ts=") ||
                boost::starts_with(key, "+ns") ||
                boost::starts_with(key, "+no_cut") ||
                boost::starts_with(key, "+h=") ||
                boost::starts_with(key, "+W=") ||
                boost::starts_with(key, "+units=") ||
                boost::starts_with(key, "+zone=")) {
      output_strings.push_back(key);
    } else if (boost::starts_with(key, "+ellps=") ||
                boost::starts_with(key, "+towgs84=") ||
                boost::starts_with(key, "+datum=")) {
      // We put these in the proj4_str for the Datum class.
      datum_strings.push_back(key);
    }
  }
  std::string strm = boost::join(output_strings, " ");
  
  // If the file contains no projection related information, we
  // supply proj.4 with a "default" interpretation that the file
  // is in geographic (unprojected) coordinates.
  if (output_strings.empty())
    m_proj_projection_str = "+proj=longlat";
  else
    m_proj_projection_str  = strm;
  
  // TODO(oalexan1): This must be done by querying the spatial ref
  if (m_proj_projection_str.find("+proj=longlat") == 0)
    m_is_projected = false;
  else
    m_is_projected = true;

  // TODO(oalexan1): Wipe this
  // TODO(oalexan1): What to do about warping?
  // Disable -180 to 180 longitude wrapping in proj4.
  // - With wrapping off, Proj4 can work significantly outside those ranges (though there is a limit)
  // - We will make sure that the input longitudes are in a safe range.
  if ((m_proj_projection_str.find("+over") == std::string::npos) &&
        (m_proj_projection_str.find("+proj=utm") == std::string::npos))
    m_proj_projection_str.append(" +over");
}

std::vector<double> GeoReference::get_towgs84_values(std::string const& s) {
  std::vector<double> o;
  std::string sub;
  if (!extract_proj4_value(s, "+towgs84", sub))
    return o;

  o.resize(6);
  int count = sscanf(sub.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf",
                     &o[0], &o[1], &o[2], &o[3], &o[4], &o[5]);
  if (count != 6)
    vw_throw(LogicErr() << "Error parsing +towgs84 from string: " << s);
  return o;
}

// Get the wkt string from the georef. It only has projection and datum information.
std::string GeoReference::get_wkt() const {
  return ogr_wkt(m_gdal_spatial_ref);
}

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
/// that point in unprojected (Geographic) coordinates (lon,lat).
Vector2 GeoReference::point_to_lonlat(Vector2 const& loc) const {
  
  if (!m_is_projected) 
    return loc;
  
  if (!m_proj_context.is_initialized())
    vw::vw_throw(vw::ArgumentErr() << "Attempted to project without a valid transform.\n");
  
  double x = loc[0];
  double y = loc[1];
  if (!m_proj_context.m_proj_to_lonlat->Transform(1, &x, &y))
    vw::vw_throw(vw::ArgumentErr() << "Failed to project point.\n");
  
  return Vector2(x, y);
}

/// Given a position in geographic coordinates (lon,lat), compute
/// the location in the projected coordinate system.
Vector2 GeoReference::lonlat_to_point(Vector2 lon_lat) const {

  if (!m_image_ll_box.empty()) {
    // Adjust lonlat to be as close as possible to the center of the image
    double mid = (m_image_ll_box.min().x() + m_image_ll_box.max().x())/2.0;
    double diff1 = std::abs(lon_lat[0] - mid);
    double diff2 = std::abs(lon_lat[0] - mid - 360);
    if (diff2 < diff1)
       lon_lat[0] -= 360;
    diff1 = std::abs(lon_lat[0] - mid);
    diff2 = std::abs(lon_lat[0] - mid + 360);
    if (diff2 < diff1)
      lon_lat[0] += 360;
  }
  
  if (!m_is_projected) 
    return lon_lat;

  if (!m_proj_context.is_initialized())
    vw::vw_throw(vw::ArgumentErr() << "Attempted to project without a valid transform.\n");
   
  double x = lon_lat[0];
  double y = lon_lat[1];
  
  // TODO(oalexan1): Must we ensure that the longitude is in the range [-180, 180]?
  
  if (!m_proj_context.m_lonlat_to_proj->Transform(1, &x, &y))
    vw::vw_throw(vw::ArgumentErr() << "Failed to project point.\n");
  return Vector2(x, y);
}

/// Convert lon/lat/alt to projected x/y/alt 
Vector3 GeoReference::geodetic_to_point(Vector3 llh) const {

  if (boost::math::isnan(llh[2]))
    return llh;
  
  Vector2 xy = lonlat_to_point(vw::math::subvector(llh, 0, 2));
  return Vector3(xy[0], xy[1], llh[2]);
}
  
/// Convert projected x/y/alt lon/lat/alt
Vector3 GeoReference::point_to_geodetic(Vector3 point) const {

  if (boost::math::isnan(point[2]))
    return point;
  
  Vector2 ll = point_to_lonlat(vw::math::subvector(point, 0, 2));
  return Vector3(ll[0], ll[1], point[2]);
}
  
// Functions for class ProjContext

// TODO(oalexan1): Wipe this!
char** ProjContext::split_proj4_string(std::string const& proj4_str, int &num_strings) {
  // TODO(oalexan1): May need to wipe this!
  std::vector<std::string> arg_strings;
  std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
  boost::split(arg_strings, trimmed_proj4_str, boost::is_any_of(" "));

  char** strings = new char*[arg_strings.size()];
  for (size_t i = 0; i < arg_strings.size(); ++i) {
    strings[i] = new char[2048];
    strncpy(strings[i], arg_strings[i].c_str(), 2048);
  }
  num_strings = boost::numeric_cast<int>(arg_strings.size());
  return strings;
}

ProjContext::ProjContext(): m_lonlat_to_proj(NULL), m_proj_to_lonlat(NULL), m_init(false) {
}

void ProjContext::init_transforms() {
  
  m_lonlat_to_proj = OGRCreateCoordinateTransformation(&m_lonlat_crs, &m_proj_crs);
  m_proj_to_lonlat = OGRCreateCoordinateTransformation(&m_proj_crs, &m_lonlat_crs);
  m_init = true;
  
  if (!m_lonlat_to_proj || !m_proj_to_lonlat)
    vw_throw(ArgumentErr() << "Failed to create the coordinate transformations "
              << "between lon-lat and projected coordinates.\n");
}

// Copy constructor. This ensures the shared pointers are recreated.
ProjContext::ProjContext(ProjContext const& other): 
  m_lonlat_crs(other.m_lonlat_crs), m_proj_crs(other.m_proj_crs), m_init(other.m_init) {

  if (!m_init) 
    return;
  
  // Recreate the transforms
  this->init_transforms();  
}

// Assignment operator
ProjContext & ProjContext::operator=(ProjContext const& other) {
  m_lonlat_crs = other.m_lonlat_crs;
  m_proj_crs = other.m_proj_crs; 
  m_init = other.m_init;

  // Recreate the transforms
  if (m_init) 
    this->init_transforms();  

  return *this;
}

      
/// Return true if the object is fully initialized
bool ProjContext::is_initialized() const {
  return m_init;
}

ProjContext::~ProjContext() {
  if (m_init) {
    // Free up the memory
    OGRCoordinateTransformation::DestroyCT(m_lonlat_to_proj);
    OGRCoordinateTransformation::DestroyCT(m_proj_to_lonlat);
  }
}
  
// Given an integer box, generate points on its boundary and the
// diagonal. We overestimate the box by ensuring the max is not exclusive.
// Otherwise, one can produce an empty box.
// TODO(oalexan1): What if the box is huge? We may need to sample it differently.
void sample_int_box(BBox2 const& pixel_bbox, std::vector<vw::Vector2> & points) {

  // Reset the output
  points.clear();

  // An empty box can have strange corners. For such, just return no points.
  if (pixel_bbox.empty()) return;
  
  // Go along the perimeter of the pixel bbox.
  for (int32 x = pixel_bbox.min().x(); x <= pixel_bbox.max().x(); x++) {
    points.push_back(Vector2(x, pixel_bbox.min().y()));
    points.push_back(Vector2(x, pixel_bbox.max().y()));
  }
  for (int32 y = pixel_bbox.min().y(); y <= pixel_bbox.max().y(); y++) {
    points.push_back(Vector2(pixel_bbox.min().x(),y));
    points.push_back(Vector2(pixel_bbox.max().x(),y));
  }
  
  // Draw an X inside the bbox. This covers the poles. It will
  // produce a lonlat boundary that is within at least one pixel of
  // the pole. This will also help catch terminator boundaries from
  // orthographic projections.
  vw::math::BresenhamLine l1(pixel_bbox.min(), pixel_bbox.max());
  while (l1.is_good()) {
    points.push_back(*l1);
    l1++;
  }
  vw::math::BresenhamLine l2(pixel_bbox.min() + Vector2i(pixel_bbox.width(), 0),
                             pixel_bbox.max() - Vector2i(pixel_bbox.width(), 0));
  while (l2.is_good()) {
    points.push_back(*l2);
    l2++;
  }
  return;
}

// Sample a float box on the edges and diagonal with a default of 100 points.
// Here the max is not assumed to be exclusive, so we sample up to and including
// the box max(). This overestimates the box a bit.
// TODO(oalexan1): 100 samples around the poles may not be good enough.
// Need to do some adaptive scheme.
void sample_float_box(BBox2 const& box, std::vector<vw::Vector2> & points,
                      int num_steps) {

  // Reset the output
  points.clear();

  // An empty box can have strange corners. For such, just return no points.
  if (box.empty()) return;

  BBox2 out_box;

  double minx = box.min().x(), maxx = box.max().x();
  double miny = box.min().y(), maxy = box.max().y();

  double rangex = maxx - minx;
  double rangey = maxy - miny;

  for (int i = 0; i <= num_steps; i++) {
    double r = double(i)/num_steps;

    // left edge
    Vector2 P2 = Vector2(minx, miny + r*rangey);
    points.push_back(P2);

    // right edge
    P2 = Vector2(maxx, miny + r*rangey);
    points.push_back(P2);

    // bottom edge
    P2 = Vector2(minx + r*rangex, miny);
    points.push_back(P2);

    // top edge
    P2 = Vector2(minx + r*rangex, maxy);
    points.push_back(P2);
    
    // diag1
    P2 = Vector2(minx + r*rangex, miny + r*rangey);
    points.push_back(P2);

    // diag2
    P2 = Vector2(maxx - r*rangex, miny + r*rangey);
    points.push_back(P2);
  }
  
}

/// For a bbox in projected space, return the corresponding bbox in pixels on
/// the image. This overestimates the box, to ensure the output box is never
/// empty if the input box is not empty.
BBox2 GeoReference::point_to_pixel_bbox(BBox2 const& point_bbox) const {
  
  // Ensure we don't get incorrect results for empty boxes with strange corners.
  if (point_bbox.empty())
  return BBox2();

  // Technically we should only have to project 2 points as the
  // georeference transform should only have a scale an translation
  // transform. Rotations are possible but outside libraries rarely
  // support it.
  BBox2 pixel_bbox;
  pixel_bbox.grow(point_to_pixel(point_bbox.min()));
  pixel_bbox.grow(point_to_pixel(point_bbox.max()));
  pixel_bbox.grow(point_to_pixel(Vector2(point_bbox.min().x(), point_bbox.max().y())));
  pixel_bbox.grow(point_to_pixel(Vector2(point_bbox.max().x(), point_bbox.min().y())));
  
  return grow_bbox_to_int(pixel_bbox);
}

// We make the max not be exclusive. This overestimates the box a bit, but
// ensures that the logic works even for subpixel boxes.
BBox2 GeoReference::pixel_to_point_bbox(BBox2 const& pixel_bbox) const {

  if (pixel_bbox.empty())
   return BBox2();
  
  BBox2 point_bbox;

  // Note that pixel_bbox().max() is exclusive.
  point_bbox.grow(pixel_to_point(pixel_bbox.min()));
  point_bbox.grow(pixel_to_point(pixel_bbox.max()));
  point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.min().x(), pixel_bbox.max().y())));
  point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.max().x(), pixel_bbox.min().y())));
  
  return point_bbox;
}

// Find the pixel box for a given lon-lat box. Since the transform from
// lon-lat to pixel is nonlinear, need to sample somewhat densely
// the edges and diagonals of the lon-lat box.
BBox2 GeoReference::pixel_to_lonlat_bbox(BBox2 const& pixel_bbox) const {

  if (pixel_bbox.empty()) return BBox2();
  
  BBox2 lonlat_bbox;
  if (!m_is_projected) {
    return pixel_to_point_bbox(pixel_bbox);
  }

  std::vector<vw::Vector2> points;
  sample_int_box(pixel_bbox, points);
  
  // Accumulate the results of sampling
  for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
    try { lonlat_bbox.grow(pixel_to_lonlat(points[ptiter])); }
    catch (const std::exception & e) {}
  }
  
  return lonlat_bbox;
}

BBox2 GeoReference::lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox, size_t nsamples) const {

  if (lonlat_bbox.empty()) return BBox2();

  if (!m_is_projected)
    return point_to_pixel_bbox(lonlat_bbox);

  BBox2 point_bbox = lonlat_to_point_bbox(lonlat_bbox, nsamples);
  return point_to_pixel_bbox(point_bbox);
}

// TODO(oalexan1): Integrate into one single function the two blocks
// having the same logic below which make use of BresenhamLine, and
// also with the function named sample_float_box().

BBox2 GeoReference::lonlat_to_point_bbox(BBox2 const& lonlat_bbox, size_t nsamples) const {
  // Alternatively this function could avoid the nsamples
  // option. The sample discrete step could just be this average
  // size of pixel in degrees.

  if (lonlat_bbox.empty()) return BBox2();

  BBox2 point_bbox;

  Vector2 lower_fraction(lonlat_bbox.width()/double(nsamples),
                          lonlat_bbox.height()/double(nsamples));
  for(size_t i = 0; i < nsamples; i++) {
      // Walk the top & bottom (technically past the edge of pixel space) rows
      double x = lonlat_bbox.min().x() + double(i) * lower_fraction.x();
      try { point_bbox.grow(lonlat_to_point(Vector2(x,lonlat_bbox.min().y()))); }
      catch (const std::exception& e) {}
      
      try { point_bbox.grow(lonlat_to_point(Vector2(x,lonlat_bbox.max().y()))); }
      catch (const std::exception& e) {}
      
      
      // Walk the left & right (technically past the edge of pixel space) columns
      double y = lonlat_bbox.min().y() + double(i) * lower_fraction.y();
      try { point_bbox.grow(lonlat_to_point(Vector2(lonlat_bbox.min().x(),y))); }
      catch (const std::exception& e) {}
      try { point_bbox.grow(lonlat_to_point(Vector2(lonlat_bbox.max().x(),y))); }
      catch (const std::exception& e) {}
  }

  // It is possible that this may not required. However in the
  // cartography it seems better to be rigorous than sorry.
  vw::math::BresenhamLine l1(Vector2i(), Vector2i(nsamples,nsamples));
  while (l1.is_good()) {
    try {
      point_bbox.grow(lonlat_to_point(elem_prod(Vector2(*l1),
                                                lower_fraction) + lonlat_bbox.min()));
    } catch (const std::exception& e) {}
    ++l1;
  }
  vw::math::BresenhamLine l2(Vector2i(nsamples,0), Vector2i(0,nsamples));
  while (l2.is_good()) {
    try {
      point_bbox.grow(lonlat_to_point(elem_prod(Vector2(*l2),lower_fraction)
                                        + lonlat_bbox.min()));
    } catch (const std::exception& e) {}
    ++l2;
  }

  return point_bbox;
}

BBox2 GeoReference::point_to_lonlat_bbox(BBox2 const& point_bbox, size_t nsamples) const {
  
  if (point_bbox.empty()) return BBox2();

  BBox2 lonlat_bbox;

  Vector2 lower_fraction(point_bbox.width()/double(nsamples),
                          point_bbox.height()/double(nsamples));

  for (size_t i = 0; i < nsamples; i++) {
    double x = point_bbox.min().x() + double(i) * lower_fraction.x();
    try { lonlat_bbox.grow(point_to_lonlat(Vector2(x,point_bbox.min().y()))); }
    catch (const std::exception& e) {}
    try { lonlat_bbox.grow(point_to_lonlat(Vector2(x,point_bbox.max().y()))); }
    catch (const std::exception& e) {}
    
    double y = point_bbox.min().y() + double(i) * lower_fraction.y();
    try { lonlat_bbox.grow(point_to_lonlat(Vector2(point_bbox.min().x(),y))); }
    catch (const std::exception& e) {}
    try { lonlat_bbox.grow(point_to_lonlat(Vector2(point_bbox.max().x(),y))); }
    catch (const std::exception& e) {}
  }
  
  // This X pattern is to capture in crossing of the poles.
  vw::math::BresenhamLine l1(Vector2i(), Vector2i(nsamples,nsamples));
  while (l1.is_good()) {
    try {
      lonlat_bbox.grow(point_to_lonlat(elem_prod(Vector2(*l1), lower_fraction)
                                          + point_bbox.min()));
    } catch (const std::exception& e) {}
    ++l1;
  }

  vw::math::BresenhamLine l2(Vector2i(nsamples,0), Vector2i(0,nsamples));
  while (l2.is_good()) {
    try {
      lonlat_bbox.grow(point_to_lonlat(elem_prod(Vector2(*l2), lower_fraction)
                                          + point_bbox.min()));
    } catch (const std::exception& e) {}
    ++l2;
  }

  return lonlat_bbox;
}

// Given a position in geographic coordinates (lon, lat), compute
// the location in pixel coordinates in this image that
// corresponds to the given geographic coordinates.
vw::Vector2 GeoReference::lonlat_to_pixel(Vector2 lon_lat) const {  
  return point_to_pixel(lonlat_to_point(lon_lat));
}

// Update the lon-lat box based on the current pixel box.
void GeoReference::ll_box_from_pix_box(BBox2 const& pixel_bbox) {
  
  if (pixel_bbox.empty())
    vw_throw(LogicErr() << "GeoReference::ll_box_from_pix_box: Empty pixel box.\n");
  
  BBox2 ll_box = pixel_to_lonlat_bbox(pixel_bbox);
  set_image_ll_box(ll_box);
  return; 
}

// The image extent lon-lat box. Used to fix 360 degree offsets.
void GeoReference::set_image_ll_box(vw::BBox2 const& bbox) {
  m_image_ll_box = bbox;
}

vw::BBox2 GeoReference::image_ll_box() const {
  return m_image_ll_box;
}

std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
  os << "-- Proj.4 Geospatial Reference Object --\n";
  if (georef.get_projcs_name() != "")
    os << "\tPROJCS name: " << georef.get_projcs_name() << "\n";
  os << "\tTransform: " << georef.transform() << "\n";
  os << "\t" << georef.datum() << "\n";
  os << "\tProj.4 String: " << georef.proj4_str() << "\n";
  
  os << "\tPixel Interpretation: ";
  if (georef.pixel_interpretation() == GeoReference::PixelAsArea)
    os << "pixel as area\n";
  else if (georef.pixel_interpretation() == GeoReference::PixelAsPoint)
    os << "pixel as point\n";
  return os;
}

}} // vw::cartography
