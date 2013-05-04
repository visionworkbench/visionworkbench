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
