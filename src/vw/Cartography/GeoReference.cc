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

#include <vw/Core/StringUtils.h>
#include <vw/Math/Geometry.h>
#include <vw/Cartography/GeoReference.h>

#if defined(VW_HAVE_PKG_GDAL)
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include "ogr_spatialref.h"
#include "cpl_string.h"
#endif

#include <vw/Math/BresenhamLine.h>
#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/FileIO/FileUtils.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

// Proj.4
#include <proj_api.h>

// Macro for checking Proj.4 output, something we do a lot of.
#define CHECK_PROJ_ERROR(ctx_input, loc) if(ctx_input.error_no()) vw_throw(ProjectionErr() << "Bad projection in GeoReference.cc. Proj.4 error: " << pj_strerrno(ctx_input.error_no()) <<".\nLocation is " << loc << ".\n")

namespace vw {
namespace cartography {

  bool read_georeference( GeoReference& georef,
                          ImageResource const& resource ) {

#if defined(VW_HAVE_PKG_GDAL)
    DiskImageResourceGDAL const* gdal =
      dynamic_cast<DiskImageResourceGDAL const*>( &resource );
    if( gdal ) 
      return read_gdal_georeference( georef, *gdal );
#endif

    DiskImageResourcePDS const* pds =
      dynamic_cast<DiskImageResourcePDS const*>( &resource );
    if( pds ) 
      return read_pds_georeference( georef, *pds );
    return false;
  }

  /// A convenience function to read georeferencing information from an image file.
  bool read_georeference( GeoReference& georef, const std::string &filename ) {

    // No image with a SPOT5 suffix can ever have georeference.
    if (vw::has_spot5_extension(filename)) return false;

    boost::shared_ptr<DiskImageResource> r(DiskImageResourcePtr( filename ));
    bool result = read_georeference( georef, *r );
    
    return result;
  }
  
  void write_georeference( ImageResource& resource,
                           GeoReference const& georef ) {
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL* gdal = dynamic_cast<DiskImageResourceGDAL*>( &resource );
    if ( gdal ) 
      return write_gdal_georeference( *gdal, georef );
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  bool read_header_string( ImageResource const& resource, std::string const& str_name,
                           std::string & str_val ) {

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL const* gdal = dynamic_cast<DiskImageResourceGDAL const*>( &resource );
    if ( gdal ) 
      return read_gdal_string( *gdal, str_name, str_val );
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  bool read_header_strings( ImageResource const& resource, 
                            std::map<std::string, std::string> & value_pairs ) {

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
    DiskImageResourceGDAL const* gdal = dynamic_cast<DiskImageResourceGDAL const*>( &resource );
    if ( gdal ) 
      return read_gdal_strings( *gdal, value_pairs );
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }

  void write_header_string( ImageResource& resource, std::string const& str_name,
                            std::string const& str_val ) {

#if defined(VW_HAVE_PKG_GDAL)
    DiskImageResourceGDAL* gdal =
      dynamic_cast<DiskImageResourceGDAL*>( &resource );
    if ( gdal ) write_gdal_string( *gdal, str_name, str_val );
    return;
#endif
    // DiskImageResourcePDS is currently read-only, so we don't bother checking for it.
    vw_throw(NoImplErr() << "This image resource does not support writing georeferencing information.");
  }


//=============================================================================================


// TODO: Move this!

size_t remove_proj4_duplicates(std::string const& str_in, std::string &str_out) {
  std::vector<std::string> arg_strings;
  std::string trimmed_proj4_str = boost::trim_copy(str_in);
  boost::split( arg_strings, trimmed_proj4_str, boost::is_any_of(" ") );
  str_out = "";
  size_t num_kept = 0;

  for (size_t i=0; i < arg_strings.size(); ++i) {
    bool duplicate = false;
    for (size_t k=i+1; k < arg_strings.size(); ++k) {
      if (arg_strings[i] == arg_strings[k]) {
        duplicate = true;
        break;
      }
    } // End k loop
    if (duplicate)
      continue;
    // Not a duplicate!
    if (num_kept > 0)
      str_out += " ";
    str_out += arg_strings[i];
    ++num_kept;
  } // End i loop

  return num_kept;
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

  void GeoReference::init_proj() {
    // Update the projection context object with the current proj4 string, 
    //  then make sure the lon center is still correct.
    m_proj_context = ProjContext( overall_proj4_str() );
    update_lon_center_private();
  }



  GeoReference::GeoReference() : m_pixel_interpretation( PixelAsArea ), m_projcs_name("") {
    set_transform(vw::math::identity_matrix<3>());
    set_geographic();
    init_proj();
  }

  GeoReference::GeoReference(Datum const& datum) :
        m_pixel_interpretation( PixelAsArea ), m_datum(datum), m_projcs_name("") {
    set_transform(vw::math::identity_matrix<3>());
    set_geographic();
    init_proj();
  }

  GeoReference::GeoReference(Datum const& datum, PixelInterpretation pixel_interpretation)
      : m_pixel_interpretation ( pixel_interpretation ), m_datum(datum), m_projcs_name("") {
    set_transform(vw::math::identity_matrix<3>());
    set_geographic();
    init_proj();
  }

  GeoReference::GeoReference(Datum const& datum,
                             Matrix<double,3,3> const& transform) :
                   m_pixel_interpretation( PixelAsArea ), m_datum(datum), m_projcs_name("") {
    set_transform(transform);
    set_geographic();
    init_proj();
  }

  GeoReference::GeoReference(Datum const& datum,
                             Matrix<double,3,3> const& transform,
                             PixelInterpretation pixel_interpretation) :
    m_pixel_interpretation(pixel_interpretation), m_datum(datum), m_projcs_name("") {
    set_transform(transform);
    set_geographic();
    init_proj();
  }

  void GeoReference::set_transform(Matrix3x3 transform) {
    m_transform = transform;
    m_shifted_transform = m_transform;
    m_shifted_transform(0,2) += 0.5*m_transform(0,0);
    m_shifted_transform(1,2) += 0.5*m_transform(1,1);
    m_inv_transform         = vw::math::inverse(m_transform);
    m_inv_shifted_transform = vw::math::inverse(m_shifted_transform);

    // If proj4 is already set up update the lon center, otherwise wait for proj4.
    if (m_proj_context.is_initialized())
      update_lon_center_private();
  }

  // We override the base classes method here so that we have the
  // opportunity to call init_proj()
  void GeoReference::set_datum(Datum const& datum) {
    m_datum = datum;

    // This is a fix for when for some reason the proj4 string
    // does not have the datum name. Example:
    // '+proj=longlat +ellps=WGS84 +no_defs '.
    if ((m_datum.spheroid_name() == "WGS_1984" ||
         m_datum.spheroid_name() == "WGS84"    ||
         m_datum.spheroid_name() == "WGS 84") &&
        (m_datum.proj4_str().find("+datum=") == std::string::npos ||
         m_datum.name() == "unknown") ){
      m_datum.name() = "WGS_1984";
      m_datum.proj4_str() += " +datum=WGS84";
    }
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

  void GeoReference::set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=stere +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +k=" << scale << " +x_0=" << false_easting
         << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_oblique_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sterea +lon_0=" << center_longitude << " +lat_0="
         << center_latitude << " +k=" << scale << " +x_0=" << false_easting
         << " +y_0=" << false_northing << " +units=m";
    set_proj4_projection_str(strm.str());
  }

  void GeoReference::set_gnomonic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=gnom +lon_0=" << center_longitude << " +lat_0="
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

    m_proj_projection_str = boost::trim_copy(s); // Store the string in this class (it is also stored in m_proj_context)

    // Extract some information from the string
    if (m_proj_projection_str.find("+proj=longlat") == 0)
      m_is_projected = false;
    else
      m_is_projected = true;

    // Disable -180 to 180 longitude wrapping in proj4.
    // - With wrapping off, Proj4 can work significantly outside those ranges (though there is a limit)
    // - We will make sure that the input longitudes are in a safe range.
    if  ( (m_proj_projection_str.find("+over") == std::string::npos) &&
          (m_proj_projection_str.find("+proj=utm") == std::string::npos) )
      m_proj_projection_str.append(" +over");

    init_proj(); // Initialize m_proj_context
    // The last step of init_proj() is to call update_lon_center_private().
  }

  void GeoReference::set_lon_center(bool centered_on_lon_zero) {
    // Don't allow switching of UTM georefs
    if (m_proj_projection_str.find("+proj=utm") != std::string::npos)
      return;
    
    // Otherwise update the field  
    m_center_lon_zero = centered_on_lon_zero;

    // Make sure that the +over flag is either in or not in the proj4
    // string as appropriate for the new lon center.
    if (m_center_lon_zero)
      clear_proj4_over();
    else 
      set_proj4_over();
  }

  bool GeoReference::safe_set_lon_center(bool new_center_around_zero) {
    bool current_center = is_lon_center_around_zero();
    if (current_center == new_center_around_zero)
      return false; // The center is already how we want it
      
    std::string proj = overall_proj4_str();
    if (proj.find("+proj=longlat") == std::string::npos)
      vw_throw( NoImplErr() << "safe_set_georef_center is only defined for longlat projections!" );
    
    // Shift the projection transform matrix to account for the new longitude center.
    Matrix3x3 affine = transform();
    if (current_center) {
      affine(0,2) += 360.0;
      set_lon_center(false);
    }
    else {
      affine(0,2) -= 360.0;
      set_lon_center(true);
    }
    set_transform(affine);
    return true; // Center was changed.
  }

  bool GeoReference::extract_proj4_value(std::string const& proj4_string, std::string const& key,
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
  bool GeoReference::extract_proj4_value(std::string const& proj4_string, std::string const& key,
                                         double &value) {
    std::string s;
    if (!extract_proj4_value(proj4_string, key, s))
      return false;
    value = atof(s.c_str());
    return true;
  }

  // Strip the "+over" text from our stored proj4 info, but don't update_lon_center().
  // - Used to strip an extra tag out of [-180,180] range images where it is not needed.
  void GeoReference::clear_proj4_over() {
    
    // Clear out m_proj_projection_str, then recreate the ProjContext object.
    if (string_replace(m_proj_projection_str, "+over", "")) {
      // If we had to make any changes, strip out any double spaces and 
      //  trailing spaces and then update our ProjContext object.
      string_replace(m_proj_projection_str, "  ", " ");
      m_proj_projection_str = boost::trim_copy(m_proj_projection_str);
      m_proj_context        = ProjContext( overall_proj4_str() );
    }
  }

  // Add the "+over" text to our stored proj4 info, but don't update_lon_center().
  void GeoReference::set_proj4_over() {
    
    // Clear out m_proj_projection_str, then recreate the ProjContext object.
    if (m_proj_projection_str.find("+over") == std::string::npos) {
      m_proj_projection_str.append(" +over");
    
      // If we had to make any changes, strip out any double spaces and 
      //  trailing spaces and then update our ProjContext object.
      string_replace(m_proj_projection_str, "  ", " ");
      m_proj_projection_str = boost::trim_copy(m_proj_projection_str);
      m_proj_context        = ProjContext( overall_proj4_str() );
    }
  }


  void GeoReference::update_lon_center(BBox2 const& pixel_bbox) {
  
    // The goal of this function is to determine which of the two standard longitude ranges
    //  ([-180 to 180] or [0 to 360]) fully contains the projected coordinate space.
    
    // UTM projections always center on 0.
    if (m_proj_projection_str.find("+proj=utm") != std::string::npos) {
      m_center_lon_zero = true;
      //std::cout << "UTM projections always center on 0.\n";
      clear_proj4_over();
      return;
    }

    // Ortho projections are tricky because pixel 0,0 may not project.
    // - Pick the longitude range where the center is closer to the projection center.
    if (m_proj_projection_str.find("+proj=ortho") != std::string::npos) {
      double lon0=0;
      m_center_lon_zero = true;
      if (extract_proj4_value(m_proj_projection_str, "+lon_0", lon0)) {
        // If the projection center is closer to 180 than it is to 0,
        //  set 180 as the projection center.
        double diff0   = math::degree_diff(lon0,   0);
        double diff180 = math::degree_diff(lon0, 180);
        if (diff180 < diff0) {
          //std::cout << "Setting ortho projection center around 180.\n";
          m_center_lon_zero = false;
          set_proj4_over();
        }
      }
      if (m_center_lon_zero)
        clear_proj4_over();
      return;
    }
    
    // Albers equal-area conic
    if (m_proj_projection_str.find("+proj=aea") != std::string::npos){
      // This seems to work best with the -180 to 180 center, so
      //  use that unless the commanded center lon is outside that range.
      double lon0 = 0;
      m_center_lon_zero = true;
      if (extract_proj4_value(m_proj_projection_str, "+lon_0", lon0)) {
        if (lon0 > 180) {
          m_center_lon_zero = false;
          set_proj4_over();
        }
      }
      if (m_center_lon_zero)
        clear_proj4_over();
      return;
    }    
    
    // See where the four corners of the image bbox project to
    std::vector<Vector2> corner_pixels;
    if (pixel_bbox.empty()) // No info, just use pixel 0,0
      corner_pixels.push_back(Vector2(0,0));
    else { 
      // BBox provided, set up all four corners.
      corner_pixels.resize(4);
      corner_pixels[0] = pixel_bbox.min();
      corner_pixels[1] = pixel_bbox.max() - Vector2(1,1);
      corner_pixels[2] = pixel_bbox.min() + Vector2(pixel_bbox.width()-1,0);
      corner_pixels[3] = pixel_bbox.min() + Vector2(0, pixel_bbox.height()-1);
    }
    
    //std::cout << "Converting pixel corners...\n";
    bool negLon = false, overLon=false;
    double minLon=99999, maxLon=-99999;
    for (size_t i=0; i<corner_pixels.size(); ++i) {
    
      // Figure out where the pixel transforms to in lon/lat.
      // - It is important that we do not normalize here!
      Vector2 point   = pixel_to_point(corner_pixels[i]);
      Vector2 lon_lat = point_to_lonlat_no_normalize(point);
      double lon = lon_lat[0]; 
      
      //printf("%d = %lf\n", i, lon);
      if (lon < minLon) minLon = lon;
      if (lon > maxLon) maxLon = lon;

      // Record the where the lonlat coordinates fall
      if (lon > 180)
        overLon = true;
      if (lon < 0)
        negLon = true;
    } // End loop through corners

    if (overLon && !negLon) { // Lons over 180, none under 0, must be 180 centered.
      m_center_lon_zero = false;
      set_proj4_over();
      return;
    }
    if (negLon && !overLon) { // Lons under 0, none over 180, must be zero centered.
      m_center_lon_zero = true;
      clear_proj4_over();
      return;
    }
    // Check for weird images with pixels that wrap-around more than 360 degrees!
    // - In this situation, determine by where the center lon is closer to.
    if (negLon && overLon) {
      double centerLon = (minLon + maxLon) / 2.0;
      double diff0     = abs(centerLon);
      double diff180   = abs(centerLon - 180);
      //printf("diff0 = %lf, diff180 = %lf\n", diff0, diff180);
      if (diff180 < diff0) { // 180 is closer, 
        m_center_lon_zero = false;
        set_proj4_over();
        return;
      } // Otherwise proceed and with the zero centered case.
    }

    // If we made it to here all pixels are in the 0-180 zone.
    // In this case, default to the more common -180 to 180 range.
    m_center_lon_zero = true;
    //std::cout << "Decreasing in shared zone, center on 0.\n";
    clear_proj4_over();      
    
    return; 
  } // End function update_lon_center

  void GeoReference::update_lon_center_private() {
    // Cal the bbox version with a zero size bbox
    update_lon_center(BBox2(0,0,0,0));
  } // End function update_lon_center_private

/*
bool GeoReference::check_projection_validity(Vector2i image_size) const {

    int NUM_PARTS = 10;
    int spacing_x = image_size[0] / (NUM_PARTS-1); // Make sure we go through most of the image
    int spacing_y = image_size[1] / (NUM_PARTS-1);

    // Loop through a grid in the image    
    for (int r=0; r<NUM_PARTS; ++r) {
      for (int c=0; c<NUM_PARTS; ++c) {
        Vector2 pixel(c*spacing_x,r*spacing_y);
        Vector2 point = pixel_to_point(pixel);
      }
    }
}
*/

double GeoReference::test_pixel_reprojection_error(Vector2 const& pixel) {
 
  Vector2 out_pixel = lonlat_to_pixel( pixel_to_lonlat(pixel) );
  Vector2 diff = out_pixel - pixel;
  double error = sqrt(diff.x()*diff.x() + diff.y()*diff.y());
  return error;
}


#if defined(VW_HAVE_PKG_GDAL)

  void GeoReference::set_wkt(std::string const& wkt) {
    const char *wkt_str = wkt.c_str();
    char **wkt_ptr = (char**)(&wkt_str);

    OGRSpatialReference gdal_spatial_ref;
    gdal_spatial_ref.importFromWkt(wkt_ptr);

    // If there is a PROJCS name, record it.
    const char * projcs = gdal_spatial_ref.GetAttrValue("PROJCS");
    if (projcs != NULL) {
      // Careful here, to avoid a segfault
      m_projcs_name = std::string(projcs);
    }

    // Create the datum. We will modify it later on.
    Datum datum;
    datum.set_datum_from_spatial_ref(gdal_spatial_ref);

    // Set the datum in the georef. Until now the georef may have been
    // completely invalid, so we need to do this step now to avoid
    // problems later on.  We'll keep on tweaking things and set the
    // datum again later one more time.
    this->set_datum(datum);

    // Read projection information out of the file
    char* proj_str_tmp;
    gdal_spatial_ref.exportToProj4(&proj_str_tmp);
    std::string proj4_str = proj_str_tmp;
    CPLFree( proj_str_tmp );

    std::vector<std::string> input_strings, output_strings, datum_strings;
    std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
    boost::split( input_strings, trimmed_proj4_str, boost::is_any_of(" ") );
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
      set_proj4_projection_str("+proj=longlat");
    else
      set_proj4_projection_str(strm);

    int utm_north = 0;
    int utm_zone = gdal_spatial_ref.GetUTMZone(&utm_north);
    if (utm_zone) set_UTM(utm_zone, utm_north);

    // Set the proj4 string for datum.
    std::string datum_proj4_ss =
        boost::trim_copy(boost::join(datum_strings, " "));
    // Add the current proj4 string in the case that our ellipse/datum
    // values are empty.
    if (datum_proj4_ss.empty()) datum_proj4_ss = datum.proj4_str();
    datum.proj4_str() = datum_proj4_ss;

    // Setting the fully processed datum
    set_datum(datum);
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

    // Create an OGRSpatialReference gdal object, load it with the
    //  proj4 string and datum information, and then use it to 
    //  generate the WKT string.

    OGRSpatialReference gdal_spatial_ref;
    Datum const& datum = this->datum();
    const std::string proj_string = this->overall_proj4_str();
    gdal_spatial_ref.importFromProj4(proj_string.c_str());

    // Apply projcs override if it was specified
    std::string projcs_name = this->get_projcs_name();
    if (!projcs_name.empty())
      gdal_spatial_ref.SetProjCS(projcs_name.c_str());

    // For perfect spheres, we set the inverse flattening to
    // zero. This is making us compliant with OpenGIS Implementation
    // Specification: CTS 12.3.10.2. In short, we are not allowed to
    // write infinity as most tools, like ArcGIS, can't read that.

    // TODO: PROJCS is still not written correctly sometimes, see
    // StereoPipelineTest/ss_mapproject_ctx_bug.

    // We also cannot handle: 
    
    // TODO: Test this some more. This is a fix for PROJCS "CH1903 / LV03"
    std::string geog_name;
    if (projcs_name != "")
      geog_name = datum.name();
    else
      geog_name = "Geographic Coordinate System";

    gdal_spatial_ref.SetGeogCS( geog_name.c_str(),
                                datum.name().c_str(),
                                datum.spheroid_name().c_str(),
                                datum.semi_major_axis(),
                                datum.semi_major_axis() == datum.semi_minor_axis() ?
                                0 : datum.inverse_flattening(),
                                datum.meridian_name().c_str(),
                                datum.meridian_offset() );

    // Make sure that this gets set properly
    std::vector<double> vals = get_towgs84_values(proj_string);
    if (vals.size() == 6)
      gdal_spatial_ref.SetTOWGS84(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]);
    
    char* wkt_str_tmp;
    gdal_spatial_ref.exportToWkt(&wkt_str_tmp);
    std::string wkt_str = wkt_str_tmp;
    OGRFree(wkt_str_tmp);

    return wkt_str;
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
  
    Vector2 lon_lat;
    if ( !m_is_projected ) {
      lon_lat = loc;
    } else {

      projXY projected;
      projLP unprojected;

      projected.u = loc[0]; // Store in proj4 object
      projected.v = loc[1];

      // Call proj4 to do the conversion and check for errors.
      unprojected = pj_inv(projected, m_proj_context.proj_ptr());
      CHECK_PROJ_ERROR(m_proj_context, loc);

      // Convert from radians to degrees.
      lon_lat = Vector2(unprojected.u * RAD_TO_DEG, unprojected.v * RAD_TO_DEG);
    }

    // Get the longitude into the correct range for this georeference.    
    lon_lat[0] = math::normalize_longitude(lon_lat[0], m_center_lon_zero);
    return lon_lat;
  }


  Vector2 GeoReference::point_to_lonlat_no_normalize(Vector2 loc) const {
    if ( !m_is_projected ) 
      return loc;

    projXY projected;
    projLP unprojected;

    projected.u = loc[0]; // Store in proj4 object
    projected.v = loc[1];

    // Call proj4 to do the conversion and check for errors.
    unprojected = pj_inv(projected, m_proj_context.proj_ptr());
    CHECK_PROJ_ERROR(m_proj_context, loc);

    // Convert from radians to degrees.
    return Vector2 (unprojected.u * RAD_TO_DEG, unprojected.v * RAD_TO_DEG);
  }

  /// Given a position in geographic coordinates (lat,lon), compute
  /// the location in the projected coordinate system.
  Vector2 GeoReference::lonlat_to_point(Vector2 lon_lat) const {

    // Get the longitude into the correct range for this georeference.    
    lon_lat[0] = math::normalize_longitude(lon_lat[0], m_center_lon_zero);

    if ( ! m_is_projected ) 
      return lon_lat;

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
    CHECK_PROJ_ERROR(m_proj_context, lon_lat);

    return Vector2(projected.u, projected.v);
  }

  /// Convert lon/lat/alt to projected x/y/alt 
  Vector3 GeoReference::geodetic_to_point(Vector3 llh) const {

    if ( boost::math::isnan(llh[2]) )
      return llh;
    
    Vector2 xy = lonlat_to_point( vw::math::subvector(llh, 0, 2) );
    return Vector3( xy[0], xy[1], llh[2] );
  }
  
  /// Convert projected x/y/alt lon/lat/alt
  Vector3 GeoReference::point_to_geodetic(Vector3 point) const {

    if ( boost::math::isnan(point[2]) )
      return point;
    
    Vector2 ll = point_to_lonlat( vw::math::subvector( point, 0, 2 ) );
    return Vector3( ll[0], ll[1], point[2] );
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


  ProjContext::ProjContext(std::string const& proj4_str ) : m_proj4_str(proj4_str) {
    m_proj_ctx_ptr.reset(pj_ctx_alloc(),pj_ctx_free);
    int num;
    char** proj_strings = split_proj4_string(m_proj4_str, num);
    m_proj_ptr.reset(pj_init_ctx( m_proj_ctx_ptr.get(),
                                  num, proj_strings ),
                     pj_free);

    VW_ASSERT( !pj_ctx_get_errno(m_proj_ctx_ptr.get()),
               InputErr() << "Proj.4 failed to initialize on string: " << m_proj4_str << "\n\tError was: " 
                          << pj_strerrno(pj_ctx_get_errno(m_proj_ctx_ptr.get())) );

    for ( int i = 0; i < num; i++ )
      delete [] proj_strings[i];
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
               InputErr() << "Proj.4 failed to initialize on string: " << m_proj4_str << "\n\tError was: " 
                          << pj_strerrno(pj_ctx_get_errno(m_proj_ctx_ptr.get())) );

    for ( int i = 0; i < num; i++ )
      delete [] proj_strings[i];
    delete [] proj_strings;
  }

  int ProjContext::error_no() const {
    return pj_ctx_get_errno(m_proj_ctx_ptr.get());
  }

//************** End functions for class ProjContext ******************
//*********************************************************************





  using vw::math::BresenhamLine;

  /// For a bbox in projected space, return the corresponding bbox in
  /// pixels on the image
  BBox2i GeoReference::point_to_pixel_bbox(BBox2 const& point_bbox) const {
    
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

  BBox2 GeoReference::pixel_to_point_bbox(BBox2i const& pixel_bbox) const {

    if (pixel_bbox.empty()) return BBox2();
    
    BBox2 point_bbox;

    // Note that pixel_bbox().max() is exclusive.
    point_bbox.grow(pixel_to_point(pixel_bbox.min()));
    point_bbox.grow(pixel_to_point(pixel_bbox.max() - Vector2(1, 1)));
    point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.min().x(), pixel_bbox.max().y()-1)));
    point_bbox.grow(pixel_to_point(Vector2(pixel_bbox.max().x()-1, pixel_bbox.min().y())));
    return point_bbox;
  }

  BBox2 GeoReference::pixel_to_lonlat_bbox(BBox2i const& pixel_bbox) const {

    if (pixel_bbox.empty()) return BBox2();
    
    BBox2 lonlat_bbox;
    if (!m_is_projected) {
      return pixel_to_point_bbox(pixel_bbox);
    }

    // Need to worry about what happens around poles
    std::vector<vw::Vector2> points;
    gen_bd_and_diag_pts(pixel_bbox, points);
    
    // Go along the perimeter of the pixel bbox.
    for (size_t ptiter = 0; ptiter < points.size(); ptiter++) {
      try { lonlat_bbox.grow(pixel_to_lonlat( points[ptiter] )); }
      catch ( const std::exception & e ) {}
    }
    
    return lonlat_bbox;
  }

  BBox2i GeoReference::lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox, size_t nsamples) const {

    if (lonlat_bbox.empty()) return BBox2();

    if (!m_is_projected) {
      return point_to_pixel_bbox(lonlat_bbox);
    }
    BBox2 point_bbox = lonlat_to_point_bbox( lonlat_bbox, nsamples );
    return point_to_pixel_bbox( point_bbox );
  }

  BBox2 GeoReference::lonlat_to_point_bbox( BBox2 const& lonlat_bbox, size_t nsamples ) const {
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
        catch ( const std::exception& e ) {}
        
        try { point_bbox.grow(lonlat_to_point(Vector2(x,lonlat_bbox.max().y()))); }
        catch ( const std::exception& e ) {}
        
        
        // Walk the left & right (technically past the edge of pixel space) columns
        double y = lonlat_bbox.min().y() + double(i) * lower_fraction.y();
        try { point_bbox.grow(lonlat_to_point(Vector2(lonlat_bbox.min().x(),y))); }
        catch ( const std::exception& e ) {}
        try { point_bbox.grow(lonlat_to_point(Vector2(lonlat_bbox.max().x(),y))); }
        catch ( const std::exception& e ) {}
    }

    // It is possible that this may not required. However in the
    // cartography it seems better to be rigorous than sorry.
    BresenhamLine l1( Vector2i(), Vector2i(nsamples,nsamples) );
    while ( l1.is_good() ) {
      try {
        point_bbox.grow( lonlat_to_point( elem_prod(Vector2(*l1),lower_fraction) + lonlat_bbox.min() ) );
      } catch ( const std::exception& e ) {}
      ++l1;
    }
    BresenhamLine l2( Vector2i(nsamples,0), Vector2i(0,nsamples) );
    while ( l2.is_good() ) {
      try {
        point_bbox.grow( lonlat_to_point( elem_prod(Vector2(*l2),lower_fraction)
                                          + lonlat_bbox.min() ) );
      } catch ( const std::exception& e ) {}
      ++l2;
    }

    return point_bbox;
  }

  BBox2 GeoReference::point_to_lonlat_bbox(BBox2 const& point_bbox, size_t nsamples) const {
    
    if (point_bbox.empty()) return BBox2();

    BBox2 lonlat_bbox;

    Vector2 lower_fraction( point_bbox.width()/double(nsamples),
                            point_bbox.height()/double(nsamples) );

    for (size_t i = 0; i < nsamples; i++ ) {
      double x = point_bbox.min().x() + double(i) * lower_fraction.x();
      try { lonlat_bbox.grow( point_to_lonlat(Vector2(x,point_bbox.min().y()))); }
      catch ( const std::exception& e ) {}
      try { lonlat_bbox.grow( point_to_lonlat(Vector2(x,point_bbox.max().y()))); }
      catch ( const std::exception& e ) {}
      
      double y = point_bbox.min().y() + double(i) * lower_fraction.y();
      try { lonlat_bbox.grow( point_to_lonlat(Vector2(point_bbox.min().x(),y))); }
      catch ( const std::exception& e ) {}
      try { lonlat_bbox.grow( point_to_lonlat(Vector2(point_bbox.max().x(),y))); }
      catch ( const std::exception& e ) {}
    }
    
    // This X pattern is to capture in crossing of the poles.
    BresenhamLine l1( Vector2i(), Vector2i(nsamples,nsamples) );
    while ( l1.is_good() ) {
      try {
        lonlat_bbox.grow( point_to_lonlat( elem_prod(Vector2(*l1), lower_fraction)
                                           + point_bbox.min() ) );
      } catch ( const std::exception& e ) {}
      ++l1;
    }

    BresenhamLine l2( Vector2i(nsamples,0), Vector2i(0,nsamples) );
    while ( l2.is_good() ) {
      try {
        lonlat_bbox.grow( point_to_lonlat( elem_prod(Vector2(*l2), lower_fraction)
                                           + point_bbox.min() ) );
      } catch ( const std::exception& e ) {}
      ++l2;
    }

    return lonlat_bbox;
  }



  std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
    os << "-- Proj.4 Geospatial Reference Object --\n";
    if (georef.get_projcs_name() != "")
      os << "\tPROJCS name: " << georef.get_projcs_name() << "\n";
    os << "\tTransform  : " << georef.transform() << "\n";
    os << "\t" << georef.datum() << "\n";
    os << "\tProj.4 String: " << georef.proj4_str() << "\n";
    os << "\tPixel Interpretation: ";
    if (georef.pixel_interpretation() == GeoReference::PixelAsArea)
      os << "pixel as area\n";
    else if (georef.pixel_interpretation() == GeoReference::PixelAsPoint)
      os << "pixel as point\n";
    if (georef.is_lon_center_around_zero())
      os << "longitude range: [-180, 180]\n";
    else
      os << "longitude range: [0, 360]\n";
    return os;
  }

  // Given an integer box, generate points on its boundary and the
  // diagonal. It is important to note that the maximum is exclusive.
  // This is used as a way of sampling the lon-lat values of all pixel values
  // in this box.
  void gen_bd_and_diag_pts(BBox2i const& pixel_bbox, std::vector<vw::Vector2> & points) {

    // Reset the output
    points.clear();

    // An empty box can have strange corners. For such, just return no points.
    if (pixel_bbox.empty()) return;
    
    // Go along the perimeter of the pixel bbox.
    for ( int32 x=pixel_bbox.min().x(); x<pixel_bbox.max().x(); ++x ) {
      points.push_back(Vector2(x,pixel_bbox.min().y()));
      points.push_back(Vector2(x,pixel_bbox.max().y()-1));
    }
    for ( int32 y=pixel_bbox.min().y()+1; y<pixel_bbox.max().y()-1; ++y ) {
      points.push_back(Vector2(pixel_bbox.min().x(),y));
      points.push_back(Vector2(pixel_bbox.max().x()-1,y));
    }
    
    // Draw an X inside the bbox. This covers the poles. It will
    // produce a lonlat boundary that is within at least one pixel of
    // the pole. This will also help catch terminator boundaries from
    // orthographic projections. Note that pixel_bbox.max() is exclusive,
    // we stop on the line right before we reach this point.
    BresenhamLine l1( pixel_bbox.min(), pixel_bbox.max() );
    while ( l1.is_good() ) {
      points.push_back(*l1);
      ++l1;
    }

    // Notice how we subtract 1 in two places to make pixel_bbox.max() exclusive.
    BresenhamLine l2( pixel_bbox.min() + Vector2i(pixel_bbox.width() - 1, 0),
                      pixel_bbox.max() - Vector2i(pixel_bbox.width(),     1) );
    while ( l2.is_good() ) {
      points.push_back(*l2 );
      ++l2;
    }

  }
  

}} // vw::cartography

#undef CHECK_PROJ_ERROR
