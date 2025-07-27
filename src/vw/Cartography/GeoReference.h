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

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/Cartography/Datum.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/Exception.h>

#include <ogr_spatialref.h>

// Boost
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

// Proj forward declarations
typedef struct pj_ctx PJ_CONTEXT;
typedef struct PJconsts PJ;  

namespace vw {
namespace cartography {

  // Define a specific exception for proj to throw.  It's derived from
  // ArgumentErr both because that's what used to be thrown here, and also
  // because basically every error proj.4 returns is due to some variety of bad input.
  VW_DEFINE_EXCEPTION(ProjectionErr, ArgumentErr);

  // This class ensures that the GeoReference transforms are properly
  // duplicated on copying.
  class ProjContext {

    /// 
    char** split_proj4_string(std::string const& proj4_str, int &num_strings);

  public:

    ProjContext();
    // The copy constructor and assignment operator are essential
    ProjContext(ProjContext const& other);
    ProjContext & operator=(ProjContext const& other);
    void init_transforms();
    ~ProjContext();

    // Transform from lonlat to projected space and their CRS
    OGRSpatialReference m_lonlat_crs, m_proj_crs;
    //boost::shared_ptr<OGRCoordinateTransformation> m_lonlat_to_proj;
    //boost::shared_ptr<OGRCoordinateTransformation> m_proj_to_lonlat;
    OGRCoordinateTransformation * m_lonlat_to_proj;
    OGRCoordinateTransformation * m_proj_to_lonlat;
    bool m_init;

    /// Return true if the object is fully initialized
    bool is_initialized() const;
  };

  // Would it make more sense for this class to keep information in an
  // OGRSpatialReference gdal object instead of in a proj4 string?
  // More information can be stored that way.

  /// The georeference class contains the mapping from image
  /// coordinates (u,v) to geospatial coordinates (typically lat/lon,
  /// or possibly meters in a UTM grid cell, etc.).  It must also
  /// encode how to translate between this coordinate system and the
  /// "Geographic" coordinate system (lon,lat)
  class GeoReference{
  
  public:
    /// The affine transform converts from pixel space to geographic
    /// or projected space and vice versa.  Most often, this process
    /// entails interpolating based on floating point pixel
    /// coordinates in the image.  However, images are discrete
    /// samples of pixel space, so you must adopt a convention
    /// regarding how floating point pixel coordinates in your
    /// georeferenced image are to be interpreted.
    ///
    /// You have one of two choices: If you assume PixelAsArea, the
    /// upper left hand corner of the top left pixel is considered as
    /// the origin (0,0), and the center of the top left pixel is
    /// (0.5, 0.5).  This assumption is common when dealing with
    /// satellite imagery or maps.
    ///
    /// On the other hand, if you assume the PixelAsPoint, then the
    /// center of the upper left hand pixel is the origin (0,0), and
    /// the top left corner of that pixel is at (-0.5,-0.5) in pixel
    /// coordinates.  This mode is common when working with elevation
    /// data, etc.
    ///
    /// Note: The Vision Workbench *always* interprets floating point
    /// pixel location (0,0) as being at the _center_ of the upper
    /// left hand pixel.  If you choose the PixelAsArea option for
    /// this flag, the GeoTransform class will automatically adjust
    /// your affine transform by (0.5,0.5) to bring the coordinate
    /// system in line with the Vision Workbench internal representation.
    ///
    /// The default pixel interpretation for GeoReference is PixelAsArea
    enum PixelInterpretation { PixelAsArea = 0, PixelAsPoint = 1 };

    // This in alias to the underlying spatial ref
    OGRSpatialReference& gdal_spatial_ref() { return m_gdal_spatial_ref; }
    OGRSpatialReference const& gdal_spatial_ref() const { return m_gdal_spatial_ref; }
   
  private:  
    PixelInterpretation m_pixel_interpretation;
    Datum       m_datum;
    Matrix<double,3,3> m_transform, m_inv_transform, m_shifted_transform, m_inv_shifted_transform;
    std::string m_proj_projection_str; // Duplicate of information in m_proj_context
    ProjContext m_proj_context;// TODO(oalexan1): Wipe this
    
    // The CRS for the georeference
    OGRSpatialReference m_gdal_spatial_ref;
    
    bool m_is_projected; // As opposed to lonlat
    
    // The image extent in projected coordinates. Used for geotransforms.
    vw::BBox2 m_image_ll_box;
    
    std::string m_projcs_name; ///< Override the projcs name when writing WKT and to file.
    
    // ---- Functions -------------------------

    /// Initialize m_proj_context with current proj4 string.
    void init_proj();
    
    void clear_proj4_over(); ///< Clears the "+over" tag from our proj4 string.
    void set_proj4_over  (); ///< Adds   the "+over" tag from our proj4 string.

    /// Attempts to extract the value of a key= part of the proj4 string.
    static bool extract_proj4_value(std::string const& proj4_string, std::string const& key,
                                    std::string &s);
    /// Attempts to extract the value of a key= part of the proj4 string and converts to double.
    static bool extract_proj4_value(std::string const& proj4_string, std::string const& key,
                                    double &value);

    /// Parse out a six element +towgs84 value, or an empty vector if it is not found.
    static std::vector<double> get_towgs84_values(std::string const& s);
    
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

    /// Takes a geodetic datum and pixel interpretation
    GeoReference(Datum const& datum, PixelInterpretation pixel_interpretation);

    /// Takes a geodetic datum and an affine transformation matrix
    GeoReference(Datum const& datum, Matrix<double,3,3> const& transform);

    /// Takes a geodetic datum and an affine transformation matrix and
    /// pixel interpretation
    GeoReference(Datum const& datum, Matrix<double,3,3> const& transform, PixelInterpretation pixel_interpretation);

    /// Destructor
    ~GeoReference();

    void set_transform(Matrix<double,3,3> transform);
    
    // Set the datum. This does not change any projection information.
    void set_datum(Datum const& datum);

    PixelInterpretation pixel_interpretation() const { return m_pixel_interpretation; }
    void set_pixel_interpretation(PixelInterpretation const& p) {
      m_pixel_interpretation = p; 
    }

    ///// Checks if an image with the given dimensions is properly contained in the projected space.
    ///// - In particular, all of the projected coordinates must fall in the selected longitude range.
    //bool check_projection_validity(Vector2i image_size) const;

    /// Returns the error in pixels for converting from a pixel to lonlat and back.
    double test_pixel_reprojection_error(Vector2 const& pixel);

    Datum const& datum() const { return m_datum; }

    // TODO: Why are there two string functions??????!!!!!!
          std::string  proj4_str() const;
    Matrix<double,3,3> transform() const { return m_transform; }

    /// Returns the proj.4 string of both the GeoReference and the datum,
    /// concatenated. This is what proj4_str() used to do.
    std::string overall_proj4_str() const;

    /// True if the georeference is using a projected coordinate
    /// system.  False if no projection is used (ie. we are only using lon, lat).
    bool is_projected() const { return m_is_projected; }

    /// Options include: WGS84, WGS72, NAD27, NAD83.  Note: you must
    /// call this routine before calling any of the routines below
    /// used to set the projection.
    void set_well_known_geogcs(std::string name);

    /// Set this georeference to use geographic coordinates (no projection)
    void set_geographic();

    /// Use an equirectangular projection
    void set_equirectangular      (double center_latitude = 0, double center_longitude = 0, double latitude_of_true_scale = 0, double false_easting = 0, double false_northing = 0);
    /// Use a sinusoidal projection
    void set_sinusoidal           (double center_longitude,    double false_easting = 0,    double false_northing = 0);
    /// Use mercator projection
    void set_mercator             (double center_latitude,     double center_longitude,     double latitude_of_true_scale = 0, double false_easting = 0, double false_northing = 0);
    /// Use transverse mercator projection
    void set_transverse_mercator  (double center_latitude,     double center_longitude,     double scale,             double false_easting = 0, double false_northing = 0);
    /// Use orthographic projection
    void set_orthographic         (double center_latitude,     double center_longitude,     double false_easting = 0, double false_northing = 0);
    /// Use stereographic projection
    void set_stereographic        (double center_latitude,     double center_longitude,     double scale,             double false_easting = 0, double false_northing = 0);
    /// Use oblique stereographic projection
    void set_oblique_stereographic(double center_latitude,     double center_longitude,     double scale,             double false_easting = 0, double false_northing = 0);
    /// Use gnomonic projection
    void set_gnomonic        (double center_latitude,     double center_longitude,     double scale,             double false_easting = 0, double false_northing = 0);
    /// Use Lambert Azimuthal projection
    void set_lambert_azimuthal    (double center_latitude,     double center_longitude,     double false_easting = 0, double false_northing = 0);
    /// Use Lambert (Conic) Conformal projection with two standard parallels
    void set_lambert_conformal    (double std_parallel_1,      double std_parallel_2,       double center_latitude,   double center_longitude,  double false_easting = 0, double false_northing = 0);
    /// Use Universal Transverse Mercator (UTM) projection
    void set_UTM(int zone, int north=true);

    /// Allows the user to explicitly specify a projection using
    /// proj.4 syntax.  The user should specify the projection only;
    /// the datum portion of the proj.4 string is still generated by the Datum object,
    /// which is very confusing.
    // Consider using instead asp::set_srs_string().
    void set_proj4_projection_str(std::string const& s);

    // Every function that modifies the georef must call this function
    void set_wkt(std::string const& wkt);

    // Get the wkt string from the georef. It only has projection and datum information.
    std::string get_wkt() const;
    
    /// Set a projcs name used in WKT and writing to disk.
    // TODO(oalexan1): May need to wipe this, or to have it get data from 
    // the underlying m_gdal_spatial_ref.
    void        set_projcs_name(std::string const& projcs_name) {m_projcs_name=projcs_name;}
    std::string get_projcs_name() const {return m_projcs_name;}
    
    //===============================================================================
    // The functions below here are for conversion between lonlat, point, and point 
    //  coordinates within the scope of a SINGLE GeoReference object.  It is 
    //  dangerous to use these functions to go between two GeoReferences because
    //  they do not take into account the potential difference in datum between 
    //  them.  To convert between coordinates of two GeoReference objects, 
    //  load them both into a GeoTransform class object (see GeoTransform.h).

    /// For a given pixel coordinate, compute the position of that
    /// pixel in this georeferenced space.
    Vector2 pixel_to_point(Vector2 pix) const;

    /// For a given location 'loc' in projected space, compute the
    /// corresponding pixel coordinates in the image.
    Vector2 point_to_pixel(Vector2 loc) const;

    /// For a point in the projected space, compute the position of
    /// that point in unprojected (Geographic) coordinates (lon,lat).
    Vector2 point_to_lonlat(Vector2 const& loc) const;

    /// Given a position in geographic coordinates (lon,lat), compute
    /// the location in the projected coordinate system.
    Vector2 lonlat_to_point(Vector2 lon_lat) const;

    /// Convert lon/lat/alt to projected x/y/alt 
    Vector3 geodetic_to_point(Vector3 llh) const;

    /// Convert projected x/y/alt lon/lat/alt
    Vector3 point_to_geodetic(Vector3 point) const;

    /// For a bbox in pixel coordinates, find what that bbox covers in lonlat
    BBox2 pixel_to_lonlat_bbox(BBox2 const& pixel_bbox) const;

    /// For a bbox in lonlat, find the bbox in pixel coordinates
    BBox2 lonlat_to_pixel_bbox(BBox2 const& lonlat_bbox, size_t nsamples = 100) const;
    
    /// For a given lonlat, provide it's point bbox
    BBox2  lonlat_to_point_bbox(BBox2 const& lonlat_bbox, size_t nsamples = 100) const;

    /// For a given bbox in point, provide the lonlat bounding box.
    BBox2  point_to_lonlat_bbox(BBox2 const& point_bbox, size_t nsamples = 100) const;
    
    /// For a given pixel coordinate, compute the position of that
    /// pixel in Geographic coordinates (lon, lat).
    Vector2 pixel_to_lonlat(Vector2 pix) const {
      return point_to_lonlat(pixel_to_point(pix));
    }

    /// Given a position in geographic coordinates (lon, lat), compute
    /// the location in pixel coordinates in this image that
    /// corresponds to the given geographic coordinates.
    Vector2 lonlat_to_pixel(Vector2 lon_lat) const;

    /// For a given pixel bbox, return the corresponding bbox in projected space
    BBox2  pixel_to_point_bbox(BBox2 const& pixel_bbox) const;

    /// For a bbox in projected space, return the corresponding bbox in
    /// pixels on the image
    BBox2 point_to_pixel_bbox(BBox2 const& point_bbox) const;
    
    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    template <class ViewT>
    BBox2 bounding_box(ImageViewBase<ViewT> const& view) const;

    /// Return the box that bounds the area represented by the
    /// geotransform for the dimensions of the given image.
    /// Note that this doesn't tell you whether the image takes the
    /// long path or the short path from the left longitude to the
    /// right longitude.
    ///
    /// The assumption here is that the projection is continuous.
    template <class ViewT>
    BBox2 lonlat_bounding_box(ImageViewBase<ViewT> const& view) const;
    
    // The image extent lon-lat box
    void set_image_ll_box(vw::BBox2 const& bbox);
    vw::BBox2 image_ll_box() const;
    // Compute based on pixel box
    void ll_box_from_pix_box(BBox2 const& pixel_bbox);
  };

  /// Format a GeoReference to a text stream
  std::ostream& operator<<(std::ostream& os, const GeoReference& georef);

  // Update a georeference based on an srs string and/or datum.
  // This function is more likely to remember the datum name than set_wkt().
  void set_srs_string(std::string srs_string, bool have_user_datum,
                      vw::cartography::Datum const& user_datum,
                      vw::cartography::GeoReference & georef);

  // Georeference I/O operations

  /// Read georeferencing information from an image resource.
  bool read_georeference(GeoReference& georef, ImageResource const& resource);

  /// A convenience function to read georeferencing information from an image file.
  bool read_georeference(GeoReference& georef, const std::string &filename);

  /// Write georeferencing information to an image resource.  You should
  /// generally call this prior to writing image data to the resource.
  void write_georeference( ImageResource& resource, GeoReference const& georef );

  /// A convenience function to write image data and its georeferencing information to a file.
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

  /// A function to read the value of a variable with given name from a geotiff
  /// file header. Return an empty string on failure.
  std::string read_header_string(std::string filename, std::string const& str_name);

  /// A function to read a string with given name from geotiff header
  bool read_header_string( ImageResource const& resource, std::string const& str_name,
                           std::string & str_val );

  /// A function to read all strings with given name from geotiff header
  bool read_header_strings( ImageResource const& resource,
                            std::map<std::string, std::string> & value_pairs);

  /// A function to write a a string with given name and value to geotiff header
  void write_header_string( ImageResource& resource, std::string const& str_name,
                            std::string const& str_val );

  // -----------------------------------------------------------
  // Template function definitions for the GeoReference class

  /// Get point coordinate bounding box of an image
  template <class ViewT>
  BBox2 GeoReference::bounding_box(ImageViewBase<ViewT> const& view) const {
    return pixel_to_point_bbox(vw::bounding_box(view.impl()));
  }

  /// Get lonlat coordinate bounding box of an image
  template <class ViewT>
  BBox2 GeoReference::lonlat_bounding_box(ImageViewBase<ViewT> const& view) const {
    return pixel_to_lonlat_bbox(vw::bounding_box(view.impl()));
  }

  // Given an integer box, generate points on its boundary and the
  // diagonal. We assume the max is not exclusive, to ensure even boxes
  // less than one pixel wide are sampled.
  void sample_int_box(BBox2 const& pixel_bbox, std::vector<Vector2> & points);

  // Sample a float box on the edges and diagonal with a default of 100 points.
  // Here the max is not assumed to be exclusive, so we sample up to and including
  // the box max(). This overestimates the box a bit.
  void sample_float_box(BBox2 const& box, std::vector<vw::Vector2> & points,
                        int num_steps = 100);
  
}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__
