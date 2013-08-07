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


#ifndef __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__
#define __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__

#include <vw/Cartography/GeoReferenceBase.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/Exception.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
#include <vw/Cartography/GeoReferenceDesc.pb.h>
#endif

/// \file GeoReference.h Class for converting between pixel and geo coordinates.

namespace vw {
namespace cartography {

  // Here is some machinery to keep track of an initialized proj.4
  // projection context using a smart pointer.
  //
  // Implementation of the methods to this class change based on what
  // Proj4 version is available. So unfortunately, the methods and
  // private variables are a mash up of what is needed for the 4.7 and
  // 4.8 versions.
  class ProjContext {
    boost::shared_ptr<void> m_proj_ctx_ptr; //  Only used for Proj4.8
    boost::shared_ptr<void> m_proj_ptr;
    std::string m_proj4_str;
    char** split_proj4_string(std::string const& proj4_str, int &num_strings);

  public:
    ProjContext() : m_proj4_str("") {};
    ProjContext(std::string const& proj4_str);
    ProjContext(ProjContext const& other ); // Only used for Proj4.8
    inline void* proj_ptr() const {
      VW_ASSERT( !m_proj4_str.empty(),
                 ArgumentErr() << "ProjContext: Projection not initialized." );
      return m_proj_ptr.get();
    }
    int error_no() const;
  };

  /// The georeference class contains the mapping from image coordinates
  /// (u,v) to geospatial coordinates (typically lat/lon, or possibly
  /// meters in a UTM grid cell, etc.)
  class GeoReference : public GeoReferenceBase {
    Matrix<double,3,3> m_transform, m_inv_transform, m_shifted_transform, m_inv_shifted_transform;
    std::string m_proj_projection_str, m_gml_str;
    ProjContext m_proj_context;
    bool m_is_projected;

    void init_proj();

    /// This method returns a version of the affine transform
    /// compatible with the VW standard notion that (0,0) is the
    /// center of the top left pixel.  If pixel_interpretation() is
    /// set to PixelAsArea, this method will adjust the affine
    /// transform my 0.5 pixels right and down.
    Matrix3x3 const& vw_native_transform() const;
    Matrix3x3 const& vw_native_inverse_transform() const;

  public:
    /// Construct a default georeference.  This georeference will use
    /// the identity matrix as the initial transformation matrix, and
    /// select the default datum (WGS84) and projection (geographic).
    GeoReference();

    /// Takes a geodetic datum.  The affine transform defaults to the
    /// identity matrix.
    GeoReference(Datum const& datum);

    /// Takes a geodetic datum and an affine transformation matrix
    GeoReference(Datum const& datum, Matrix<double,3,3> const& transform);

    /// Takes a geodetic datum and an affine transformation matrix and
    /// pixel interpretation
    GeoReference(Datum const& datum, Matrix<double,3,3> const& transform, PixelInterpretation pixel_interpretation);

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
    /// Construct a GeoReference from a GeoReferenceDesc
    GeoReference(GeoReferenceDesc const& desc);

    /// Create a GeoReferenceDesc for the georef
    GeoReferenceDesc build_desc();
#endif

    /// Destructor.
    virtual ~GeoReference() {}

    void set_transform(Matrix<double,3,3> transform);
    virtual void set_datum(Datum const& datum);

    std::string proj4_str() const;
    const std::string gml_str()    const { return m_gml_str; }
    Matrix<double,3,3> transform() const { return m_transform; }

    // Returns the proj.4 string of both the GeoReference and the datum,
    // concatenated. This is what proj4_str() used to do.
    std::string overall_proj4_str() const;

    /// True if the georeference is using a projected coordinate
    /// system.  False if no projection is used (ie. we are only using
    /// lon, lat).
    bool is_projected() const { return m_is_projected; }

    /// Options include: WGS84, WGS72, NAD27, NAD83.  Note: you must
    /// call this routine before calling any of the routines below
    /// used to set the projection.
    void set_well_known_geogcs(std::string name);

    /// Set this georeference to use geographic coordinates (no projection)
    void set_geographic();

    /// Use an equirectangular projection
    void set_equirectangular(double center_latitude = 0, double center_longitude = 0, double latitude_of_true_scale = 0, double false_easting = 0, double false_northing = 0);
    /// Use a sinusoidal projection
    void set_sinusoidal(double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use mercator projection
    void set_mercator(double center_latitude, double center_longitude, double latitude_of_true_scale = 0, double false_easting = 0, double false_northing = 0);
    /// Use transverse mercator projection
    void set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use orthographic projection
    void set_orthographic(double center_latitude, double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use steregraphic projection
    void set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use oblique steregraphic projection
    void set_oblique_stereographic(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use Lambert Azimuthal projection
    void set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use Lambert (Conic) Conformal projection with two standard parallels
    void set_lambert_conformal(double std_parallel_1, double std_parallel_2, double center_latitude, double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use Universal Transverse Mercator (UTM) projection
    void set_UTM(int zone, int north=true);

    /// Allows the user to explicitly specify a projection using
    /// proj.4 syntax.  The user should specify the projection only;
    /// the datum portion of the proj.4 string is still generated by
    /// the Datum object.
    void set_proj4_projection_str(std::string const& s);

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL
    // Loads the datum and projection information from the given
    // string in WKT ("Well-Known Text") format.
    void set_wkt(std::string const& wkt);
#endif

    /// For a given pixel coordinate, compute the position of that
    /// pixel in this georeferenced space.
    virtual Vector2 pixel_to_point(Vector2 pix) const;

    /// For a given location 'loc' in projected space, compute the
    /// corresponding pixel coordinates in the image.
    virtual Vector2 point_to_pixel(Vector2 loc) const;

    /// For a point in the projected space, compute the position of
    /// that point in unprojected (Geographic) coordinates (lat,lon).
    virtual Vector2 point_to_lonlat(Vector2 loc) const;

    /// Given a position in geographic coordinates (lat,lon), compute
    /// the location in the projected coordinate system.
    virtual Vector2 lonlat_to_point(Vector2 lon_lat) const;

    /// For a bbox in pixel coordinates, find what that bbox covers
    /// in lonlat
    virtual BBox2 pixel_to_lonlat_bbox(BBox2i const& pixel_bbox) const {
      if (!m_is_projected) {
        return pixel_to_point_bbox(pixel_bbox);
      }
      return GeoReferenceBase::pixel_to_lonlat_bbox(pixel_bbox);
    }

    /// For a bbox in lonlat, find the bbox in pixel coordinates
    virtual BBox2i lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox, size_t nsamples = 100) const {
      if (!m_is_projected) {
        return point_to_pixel_bbox(lonlat_bbox);
      }
      return GeoReferenceBase::lonlat_to_pixel_bbox(lonlat_bbox, nsamples);
    }
  };

  inline std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
    os << "-- Proj.4 Geospatial Reference Object --\n";
    os << "\tTransform  : " << georef.transform() << "\n";
    os << "\t" << georef.datum() << "\n";
    os << "\tProj.4 String: " << georef.proj4_str() << "\n";
    os << "\tPixel Interpretation: ";
    if (georef.pixel_interpretation() == GeoReference::PixelAsArea)
      os << "pixel as area\n";
    else if (georef.pixel_interpretation() == GeoReference::PixelAsPoint)
      os << "pixel as point\n";
    return os;
  }

  //
  // Georeference I/O operations
  //

  /// Read georeferencing information from an image resource.
  bool read_georeference( GeoReference& georef, ImageResource const& resource );

  /// A convenience function to read georeferencing information from an image file.
  inline bool read_georeference( GeoReference& georef, const std::string &filename ) {
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::open( filename ));
    bool result = read_georeference( georef, *r );
    return result;
  }

  /// A convenience function to read an image and its georeferencing information.
  template <class PixelT>
  bool read_georeferenced_image( ImageView<PixelT>& image,
                                 GeoReference& georef,
                                 const std::string &filename ) {
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::open( filename ));
    bool result = read_georeference( georef, *r );
    read_image( image, *r );
    return result;
  }

  /// Write georeferencing information to an image resource.  You should
  /// generally call this prior to writing image data to the resource.
  void write_georeference( ImageResource& resource, GeoReference const& georef );

  /// A convenience function to write image data and its georeferencing information
  /// to a file.
  template <class ImageT>
  void write_georeferenced_image( std::string const& filename,
                                  ImageViewBase<ImageT> const& image,
                                  GeoReference const& georef,
                                  ProgressCallback const& progress_callback = ProgressCallback::dummy_instance() ) {
    vw_out(InfoMessage, "fileio") << "\tSaving image: " << filename << "\t";
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create( filename, image.format() ));
    vw_out(InfoMessage, "fileio") << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";
    write_georeference( *r, georef );
    write_image( *r, image, progress_callback );
  }

  /// A function to read a string with given name from geotiff header
  bool read_header_string( ImageResource const& resource, std::string const& str_name,
                           std::string & str_val );

  /// A function to write a a string with given name and value to geotiff header
  void write_header_string( ImageResource& resource, std::string const& str_name,
                            std::string const& str_val );

  /// The following namespace contains functions that return GeoReferences
  /// for certain well-known output styles, such as KML (and related
  /// functions involved in doing so).
  namespace output {
    namespace kml {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        Vector2 pos = tx.forward( pixel );
        Vector2 x_vector = tx.forward( pixel+Vector2(1,0) ) - pos;
        Vector2 y_vector = tx.forward( pixel+Vector2(0,1) ) - pos;
        double degrees_per_pixel = (std::min)( norm_2(x_vector), norm_2(y_vector) );
        double pixels_per_circumference = 360.0 / degrees_per_pixel;
        int scale_exponent = (int) ceil( log(pixels_per_circumference)/log(2.0) );
        if (scale_exponent >= 31) scale_exponent = 30;
        return 1 << scale_exponent;
      }
    } // namespace: vw::cartography::output::kml

    namespace tms {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        // It's exactly the same as the one for KML.
        return vw::cartography::output::kml::compute_resolution(tx, pixel);
      }
    } // namespace vw::cartography::output::tms
  } // namespace vw::cartography::output

  // Simple GeoReference modification tools
  GeoReference crop( GeoReference const& input,
                     double upper_left_x, double upper_left_y,
                     double width=0, double height=0 );
  GeoReference crop( GeoReference const& input,
                     BBox2 const& bbox );

  GeoReference resample( GeoReference const& input,
                         double scale_x, double scale_y );
  GeoReference resample( GeoReference const& input,
                         double scale );

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__
