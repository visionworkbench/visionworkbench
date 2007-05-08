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
#ifndef __VW_CARTOGRAPHY_GEOREFERENCE_H__
#define __VW_CARTOGRAPHY_GEOREFERENCE_H__

#include <vw/Image/PerPixelViews.h>
#include <vw/Core/Functors.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/Datum.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/FileMetadata.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/smart_ptr.hpp>
 
namespace vw {
namespace cartography {
  
  /// The georeference class contains the mapping from image coordinates
  /// (u,v) to geospatial coordinates (typically lat/lon, or possibly
  /// meters in a UTM grid cell, etc.)
  class GeoReference : public FileMetadata {
    std::string m_name;
    Matrix<double,3,3> m_transform;
    std::string m_proj4_str, m_wkt_str, m_gml_str;
    bool m_is_projected;
    int m_pixel_interpretation;

  public:

    /// Construct a default georeference.  This georeference will use
    /// the identity matrix as the initial transformation matrix, and
    /// select the default datum (WGS84) and projection (geographic).
    GeoReference();
    
    /// Takes a string in proj.4 format. The affine transform defaults to the identity matrix.
    GeoReference(std::string const proj4_str);
    /// Takes a string in proj.4 format and an affine transformation matrix.
    GeoReference(std::string const proj4_str, Matrix<double,3,3> const& transform);

    /// Takes a void pointer to an OGRSpatialReference. The affine transform defaults to the identity matrix.
    GeoReference(void* spatial_ref_ptr);
    /// Takes a void pointer to an OGRSpatialReference and an affine transformation matrix.
    GeoReference(void* spatial_ref_ptr, Matrix<double,3,3> const& transform); 

    /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
    GeoReference(GeoDatum const& datum);
    /// Takes a geodetic datum and an affine transformation matrix
    GeoReference(GeoDatum const& datum, Matrix<double,3,3> const& transform);
    
    /// Destructor.
    virtual ~GeoReference() {}
    
    /// Implementation of FileMetadata interface.
    static std::string metadata_type_static(void) { return "GeoReference"; };
    virtual std::string metadata_type(void) const { return metadata_type_static(); };
    virtual void read_file_metadata(DiskImageResource* r);
    virtual void write_file_metadata(DiskImageResource* r) const;
    static void register_disk_image_resource(std::string const& disk_image_resource_type,
                                             read_metadata_func read_func,
                                             write_metadata_func write_func);

    /// Takes a void pointer to an OGRSpatialReference
    void set_spatial_ref(void* spatial_ref_ptr);
    void set_transform(Matrix<double,3,3> transform) { m_transform = transform; }
    void set_proj4_str(std::string const& proj4_str);
    void set_wkt_str(std::string const& wkt_str);
    
    const std::string   proj4_str()  const { return m_proj4_str; }
    const std::string   wkt_str()    const { return m_wkt_str; }
    const std::string   gml_str()    const { return m_gml_str; }
    const void*         spatial_ref_ptr() const;
    GeoDatum datum() const;
    std::string projection_name() const;
    Matrix<double,3,3> transform() const { return m_transform; }
    bool is_projected() const { return m_is_projected; }



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
    /// your affine transform my (0.5,0.5) to bring the coordinate
    /// system in line with the Vision Workbench internal
    /// representation.
    ///
    /// The default pixel interpretation for GeoReference is PixelAsPoint
    enum PixelInterpretation { PixelAsArea, PixelAsPoint };

    PixelInterpretation pixel_interpretation() const { return (PixelInterpretation)(m_pixel_interpretation); }
    void set_pixel_interpretation(PixelInterpretation const& p) { m_pixel_interpretation = p; }

    /// This method returns a version of the affine transform
    /// compatible with the VW standard notion that (0,0) is the
    /// center of the top left pixel.  If pixel_interpretation() is
    /// set to PixelAsArea, this method will adjust the affine
    /// transform my 0.5 pixels right and down.
    Matrix<double,3,3> vw_native_transform() const;

    /// Return the box that bounds the area represented by the
    /// geotransform for an image of the given dimensions.
    BBox2 bounding_box(int width, int height) const;    

    /// Options include: WGS84, WGS72, NAD27, NAD83, or EPSG:n where n
    /// is the four digit EPSG code number.  Note: you must call this
    /// routine before calling any of the routines below used to set
    /// the projection.
    void set_well_known_geogcs(std::string name);

    /// Set this georeference to use a sinusoidal projection
    void set_sinusoidal(double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use mercator projection
    void set_mercator(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use transverse mercator projection
    void set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use orthographic projection
    void set_orthographic(double center_latitude, double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use steregraphic projection
    void set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use steregraphic projection
    void set_polar_stereographic(double center_latitude, double center_longitude, double scale, double false_easting = 0, double false_northing = 0);
    /// Use Lambert Azimuthal projection
    void set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting = 0, double false_northing = 0);
    /// Use Universal Transverse Mercator (UTM) projection
    void set_UTM(int zone, int north=true);
  };
  
  inline std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
    os << "-- Geospatial Reference Object --\n";
    os << "\tTransform  : " << georef.transform() << "\n";
    os << "\t" << georef.datum() << "\n";
    os << "\tProjection: " << georef.projection_name() << "\n";
    os << "\tProj.4 String: " << georef.proj4_str() << "\n";
    return os;
  }

  template <class ElemT>
  class XYZtoLonLatFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    XYZtoLonLatFunctor(bool east_positive = true) : m_east_positive(east_positive) {}

    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
      // Deal with "missing pixels"
      if (p == Vector<ElemT,3>()) { return p; }

      double radius = norm_2(p);
      double sin_lat = p.z() / radius;
      
      double cos_lat = sqrt(1.0 - sin_lat * sin_lat);
      double lat = asin(sin_lat);
      double lon;
      if (m_east_positive) 
        lon = atan2(p.y(), p.x());
      else // West positive longitude
        lon = atan2(-p.y(), p.x()); 

      // For consistency-sake, we always return a positive longitude.
      if (lon < 0) 
        lon = 2*M_PI + lon;

      return Vector<ElemT,3> (lon * 180.0 / M_PI, lat * 180.0 / M_PI, radius);
    }
  };
  

  /// Takes an ImageView of Vector<ElemT,3> in cartesian 3 space and
  /// returns a ImageView of vectors that contains the lat, lon, and
  /// radius of that point.  For consistency with cartographic
  /// convention, angular values are return in degrees rather than
  /// radians.
  ///
  /// Note: The following assumes latitude is measured from the
  /// equatorial plane with north positive. This is different than
  /// normal spherical coordinate conversion where the equivalent
  /// angle is measured from the positive z axis.
  ///
  /// Note: notice that the order of the returned triple is longitude,
  /// latitude, radius.  This ordering of lon/lat is consistent with
  /// the notion of horizontal (x) and vertical (y) coordinates in an
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT, XYZtoLonLatFunctor<typename ImageT::pixel_type::value_type> >
  inline xyz_to_lon_lat_radius( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,XYZtoLonLatFunctor<vector_value_type> >( image.impl(), XYZtoLonLatFunctor<vector_value_type>(east_positive) );
  }


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_H__
