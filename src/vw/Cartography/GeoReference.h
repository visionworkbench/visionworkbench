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
#ifndef __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__
#define __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__

#include <vw/Core/Functors.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/GeoReferenceBase.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/FileMetadata.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/smart_ptr.hpp>
 
namespace vw {
namespace cartography {

  // Forward declaration
  class ProjContext;

  
  /// The georeference class contains the mapping from image coordinates
  /// (u,v) to geospatial coordinates (typically lat/lon, or possibly
  /// meters in a UTM grid cell, etc.)
  class GeoReference : public GeoReferenceBase, public FileMetadata {
    Matrix<double,3,3> m_transform, m_inv_transform, m_shifted_transform, m_inv_shifted_transform;
    std::string m_proj_projection_str, m_gml_str;
    boost::shared_ptr<ProjContext> m_proj_context;
    bool m_is_projected;
    
    void init_proj();

    /// This method returns a version of the affine transform
    /// compatible with the VW standard notion that (0,0) is the
    /// center of the top left pixel.  If pixel_interpretation() is
    /// set to PixelAsArea, this method will adjust the affine
    /// transform my 0.5 pixels right and down.
    Matrix<double,3,3> vw_native_transform() const;
    Matrix<double,3,3> vw_native_inverse_transform() const;

  public:

    /// Construct a default georeference.  This georeference will use
    /// the identity matrix as the initial transformation matrix, and
    /// select the default datum (WGS84) and projection (geographic).
    GeoReference();

    /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
    GeoReference(Datum const& datum);
    
    /// Takes a geodetic datum and an affine transformation matrix
    GeoReference(Datum const& datum, Matrix<double,3,3> const& transform);
    
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

    void set_transform(Matrix<double,3,3> transform);
    virtual void set_datum(Datum const& datum);

    std::string proj4_str() const;
    const std::string gml_str()    const { return m_gml_str; }
    Matrix<double,3,3> transform() const { return m_transform; }

    /// Options include: WGS84, WGS72, NAD27, NAD83.  Note: you must
    /// call this routine before calling any of the routines below
    /// used to set the projection.
    void set_well_known_geogcs(std::string name);

    /// Set this georeference to use geographic coordinates (no projection)
    void set_geographic();

    /// Set this georeference to use a sinusoidal projection
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
    /// Use Universal Transverse Mercator (UTM) projection
    void set_UTM(int zone, int north=true);
    
    /// Allows the user to explicitly specify a projection using
    /// proj.4 syntax.  The user should specify the projection only;
    /// the datum portion of the proj.4 string is still generated by
    /// the Datum object.
    void set_proj4_projection_str(std::string const& s);

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
  };
  
  inline std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
    os << "-- Proj.4 Geospatial Reference Object --\n";
    os << "\tTransform  : " << georef.transform() << "\n";
    os << "\t" << georef.datum() << "\n";
    os << "\tProj.4 String: " << georef.proj4_str() << "\n";
    return os;
  }


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_PROJGEOREFERENCE_H__
