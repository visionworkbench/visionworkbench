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
#include <vw/Cartography/Projection.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/smart_ptr.hpp>
 
namespace vw {
namespace cartography {
  
  /// The georeference class contains the mapping from image coordinates
  /// (u,v) to geospatial coordinates (typically lat/lon, or possibly
  /// meters in a UTM grid cell, etc.)
  class GeoReference {
    std::string m_name;
    Matrix<double,3,3> m_transform;
    std::string m_proj4_str, m_wkt_str;
    bool m_is_projected;

  public:
    /// Construct a default georeference.  This georeference will use
    /// the identity matrix as the initial transformation matrix, and
    /// select the default datum (???) and projection (geographic).
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

    /// Takes a void pointer to an OGRSpatialReference
    void set_spatial_ref(void* spatial_ref_ptr);
    void set_transform(Matrix<double,3,3>& transform) { m_transform = transform; }
    void set_proj4_str(std::string const& proj4_str);
    void set_wkt_str(std::string const& wkt_str);
    
    const std::string   proj4_str()  const { return m_proj4_str; }
    const std::string   wkt_str()    const { return m_wkt_str; }
    const void*         spatial_ref_ptr() const;
    GeoDatum datum() const;
    GeoProjection projection() const;
    Matrix<double,3,3> transform() const {
      return m_transform;
    }
    bool is_projected() const;

    /// Options include: WGS84, WGS72, NAD27, NAD83, or EPSG:n where n
    /// is the four digit EPSG code number.  Note: you must call this
    /// routine before calling any of the routines below used to set
    /// the projection.
    void set_well_known_geogcs(std::string name);

    /// Set this georeference to use a sinusoidal projection
    void set_sinusoidal(double center_longitude, 
                        double false_easting = 0,
                        double false_northing = 0);
  };
  
  inline std::ostream& operator<<(std::ostream& os, const GeoReference& georef) {
    os << "-- Geospatial Reference Object --\n";
    os << "\tTransform  : " << georef.transform() << "\n";
    os << "\t" << georef.datum() << "\n";
    os << "\t" << georef.projection() << "\n";
    os << "\tProj.4: " << georef.proj4_str() << "\n";
    return os;
  }


  template <class ElemT>
  class XYZtoLatLonFunctor : public UnaryReturnSameType {
    bool m_east_positive;
  public:
    XYZtoLatLonFunctor(bool east_positive = true) : m_east_positive(east_positive) {}

    Vector<ElemT,3> operator()(Vector<ElemT,3> const& p) const {
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

      return Vector<ElemT,3> (lat * 180.0 / M_PI, lon * 180.0 / M_PI, radius);
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
  template <class ImageT>
  UnaryPerPixelView<ImageT, XYZtoLatLonFunctor<typename ImageT::pixel_type::value_type> >
  inline xyz_to_latlon( ImageViewBase<ImageT> const& image, bool east_positive = true ) {
    typedef typename ImageT::pixel_type::value_type vector_value_type;
    return UnaryPerPixelView<ImageT,XYZtoLatLonFunctor<vector_value_type> >( image.impl(), XYZtoLatLonFunctor<vector_value_type>(east_positive) );
  }


}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_H__
