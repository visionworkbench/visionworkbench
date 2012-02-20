// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Cartography/PointImageManipulation.h>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace vw;

Vector3 cartography::GeodeticToCartesian::operator()( Vector3 const& v ) const {
  if ( boost::math::isnan(v[2]) )
    return Vector3();
  return m_datum.geodetic_to_cartesian(v);
}

Vector3 cartography::CartesianToGeodetic::operator()( Vector3 const& v ) const {
  if ( v == Vector3() )
    return Vector3(0,0,std::numeric_limits<double>::quiet_NaN());
  return m_datum.cartesian_to_geodetic(v);
}

Vector3 cartography::GeodeticToProjection::operator()( Vector3 const& v ) const {
  if ( boost::math::isnan(v[2]) )
    return v;
  Vector2 pix = m_reference.lonlat_to_pixel( subvector(v, 0, 2) );
  return Vector3( pix[0], pix[1], v[2] );
}

Vector3 cartography::ProjectionToGeodetic::operator()( Vector3 const& v ) const {
  Vector2 ll = m_reference.pixel_to_lonlat( subvector( v, 0, 2 ) );
  return Vector3( ll[0], ll[1], v[2] );
}

Vector3 cartography::GeodeticToPoint::operator()( Vector3 const& v ) const {
  if ( boost::math::isnan(v[2]) )
    return v;
  Vector2 pix = m_reference.lonlat_to_point( subvector(v, 0, 2) );
  return Vector3( pix[0], pix[1], v[2] );
}

Vector3 cartography::PointToGeodetic::operator()( Vector3 const& v ) const {
  Vector2 ll = m_reference.point_to_lonlat( subvector( v, 0, 2 ) );
  return Vector3( ll[0], ll[1], v[2] );
}
