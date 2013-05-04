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


/// \file Sphere.h
///
/// A sphere shape class.
#ifndef __VW_GEOMETRY_SPHERE_H__
#define __VW_GEOMETRY_SPHERE_H__

#include <iostream>
#include <cmath>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>
#include <vw/Geometry/Shape.h>

namespace vw {
namespace geometry {

  // *******************************************************************
  // class SphereBase
  // *******************************************************************

  /// A general n-dimensional sphere class,
  /// represented by a vector pointing to the center, and a radius.
  template <class SphereT, class RealT, int DimN>
  class SphereBase : public ShapeBase<SphereBase<SphereT, RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty sphere.
    SphereBase() : m_empty( true ) {}

    /// Constructs a sphere with the given center and radius.
    template <class VectorT>
    SphereBase( VectorBase<VectorT> const& center, RealT radius ) :
      m_center( center ), m_radius( radius ), m_empty( false ) {}

    /// Constructs a 2D sphere with the given center point
    /// coordinates and radius.  (Only valid for 2D spheres.)
    SphereBase( RealT centerx, RealT centery, RealT radius )
      : m_center( Vector<RealT,2>(centerx,centery) ),
        m_radius( radius ), m_empty( false ) {}

    /// Constructs a 3D sphere with the given center point
    /// coordinates and radius.  (Only valid for 3D spheres.)
    SphereBase( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : m_center( Vector<RealT,3>(centerx,centery,centerz) ),
        m_radius( radius ), m_empty( false ) {}

    /// Returns the derived implementation type.
    SphereT& impl() { return *static_cast<SphereT*>(this); }

    /// Returns the derived implementation type.
    SphereT const& impl() const { return *static_cast<SphereT const*>(this); }

    /// Grows a sphere to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(m_center.size() == 0 || point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      if (empty()) {
        m_center = point;
        m_radius = 0;
        m_empty = false;
      }
      else {
        Vector<RealT, DimN> v = point - m_center;
        double dist = norm_2(v);
        double radius_d;
        if (dist > (double)m_radius) {
          radius_d = (dist - (double)m_radius)/2;
          m_center += v/dist*radius_d;
          v = point - m_center;
          dist = norm_2(v);
          if (dist > 0)
            m_radius += (RealT)dist;
        }
      }
    }

    /// Grows a sphere to include the given sphere.
    template <class SphereT1, class RealT1, int DimN1>
    void grow( SphereBase<SphereT1, RealT1, DimN1> const& sphere ) {
      if (sphere.empty())
        return;
      if (empty()) {
        m_center = sphere.center_();
        m_radius = sphere.radius();
        m_empty = sphere.empty();
      }
      else {
        Vector<RealT, DimN> v = sphere.center_() - m_center;
        grow(sphere.center_() + normalize(v)*sphere.radius());
      }
    }

    /// Crops (intersects) this sphere to the given sphere.
    template <class SphereT1, class RealT1, int DimN1>
    void crop( SphereBase<SphereT1, RealT1, DimN1> const& sphere ) {
      if (empty() || sphere.empty())
        return;
      Vector<RealT, DimN> v = sphere.center_() - m_center;
      RealT dist = (RealT)norm_2(v);
      if (dist > (m_radius + sphere.radius())) {
        // intersects(sphere) == false
        m_radius = 0;
        m_empty = true;
      }
      else if (dist <= (m_radius - sphere.radius())) {
        // contains(sphere) == true
        m_center = sphere.center_();
        m_radius = sphere.radius();
      }
      else if (dist <= (sphere.radius() - m_radius)) {
        // sphere.contains(*this) == true
        // smallest Sphere is *this
        return;
      }
      else if (m_radius >= dist) {
        // smallest Sphere is sphere
        m_center = sphere.center_();
        m_radius = sphere.radius();
      }
      else if (sphere.radius() >= dist) {
        // smallest Sphere is *this
        return;
      }
      else {
        // normal intersection
        RealT x = (dist + (m_radius*m_radius - sphere.radius()*sphere.radius())/dist)/2;
        m_center += v/dist*x;
        m_radius = std::sqrt(m_radius*m_radius - x*x);
      }
    }

    /// Expands this sphere by the given offset in every direction.
    void expand( RealT offset ) {
      if (!empty())
        m_radius += offset;
    }

    /// Contracts this sphere by the given offset in every direction.
    void contract( RealT offset ) {
      if (!empty())
        m_radius -= offset;
    }

    /// Returns true if the given point is contained in the sphere.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      VW_ASSERT(m_center.size() == 0 || point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      return (!empty() && (norm_2(m_center - point) <= (double)m_radius));
    }

    /// Returns true if the given sphere is entirely contained
    /// in this sphere.
    template <class SphereT1, class RealT1, int DimN1>
    bool contains( const SphereBase<SphereT1, RealT1, DimN1> &sphere ) const {
      return (!empty() && !sphere.empty() && (norm_2(m_center - sphere.center_()) <= (double)(m_radius - sphere.radius())));
    }

    /// Returns true if the given sphere intersects this
    /// sphere.
    template <class SphereT1, class RealT1, int DimN1>
    bool intersects( const SphereBase<SphereT1, RealT1, DimN1>& sphere ) const {
      return (!empty() && !sphere.empty() && (norm_2(m_center - sphere.center_()) <= (double)(m_radius + sphere.radius())));
    }

    /// Returns the size (i.e. the diameter) of the sphere.
    RealT size() const { return 2*m_radius; }

    /// Returns the center point of the sphere.
    Vector<RealT, DimN> center() const { return m_center; }

    /// Returns the center point of the sphere.
    Vector<RealT, DimN> const& center_() const { return m_center; }

    /// Returns the center point of the sphere.
    Vector<RealT, DimN> &center_() { return m_center; }

    /// Returns the radius of the sphere.
    RealT const& radius() const { return m_radius; }

    /// Returns the radius of the sphere.
    RealT &radius() { return m_radius; }

    /// Returns true if the sphere is empty (i.e. degenerate).
    /// Note that even a zero-radius sphere is not empty,
    /// because spheres contain their surfaces: a zer-radius
    /// sphere actually consists of a single point.
    bool empty() const { return m_empty; }

    /// Scales the sphere relative to the origin.
    template <class ScalarT>
    SphereBase& operator*=( ScalarT s ) {
      m_center *= s;
      m_radius *= s;
      return *this;
    }

    /// Scales the sphere relative to the origin.
    template <class ScalarT>
    SphereBase& operator/=( ScalarT s ) {
      m_center /= s;
      m_radius /= s;
      return *this;
    }

    /// Offsets the sphere by the given vector.
    template <class VectorT>
    SphereBase& operator+=( VectorBase<VectorT> const& v ) {
      m_center += v;
      return *this;
    }

    /// Offsets the sphere by the negation of the given vector.
    template <class VectorT>
    SphereBase& operator-=( VectorBase<VectorT> const& v ) {
      m_center -= v;
      return *this;
    }

  protected:
    Vector<RealT, DimN> m_center;
    RealT m_radius;
    bool m_empty;
  };

  /// Scales a sphere relative to the origin.
  template <class SphereT, class RealT, int DimN, class ScalarT>
  inline SphereT operator*( SphereBase<SphereT, RealT, DimN> const& sphere, ScalarT s ) {
    SphereT result = sphere.impl();
    result *= s;
    return result;
  }

  /// Scales a sphere relative to the origin.
  template <class SphereT, class RealT, int DimN, class ScalarT>
  inline SphereT operator/( SphereBase<SphereT, RealT, DimN> const& sphere, ScalarT s ) {
    SphereT result = sphere.impl();
    result /= s;
    return result;
  }

  /// Scales a sphere relative to the origin.
  template <class SphereT, class RealT, int DimN, class ScalarT>
  inline SphereT operator*( ScalarT s, SphereBase<SphereT, RealT, DimN> const& sphere ) {
    return sphere * s;
  }

  /// Offsets a sphere by the given vector.
  template <class SphereT, class RealT, int DimN, class VectorT>
  inline SphereT operator+( SphereBase<SphereT, RealT, DimN> const& sphere, VectorBase<VectorT> const& v ) {
    SphereT result = sphere.impl();
    result += v.impl();
    return result;
  }

  /// Offsets a sphere by the given vector.
  template <class SphereT, class RealT, int DimN, class VectorT>
  inline SphereT operator+( VectorBase<VectorT> const& v, SphereBase<SphereT, RealT, DimN> const& sphere ) {
    return sphere + v;
  }

  /// Offsets a sphere by the negation of the given vector.
  template <class SphereT, class RealT, int DimN, class VectorT>
  inline SphereT operator-( SphereBase<SphereT, RealT, DimN> const& sphere, VectorBase<VectorT> const& v ) {
    SphereT result = sphere.impl();
    result -= v.impl();
    return result;
  }

  /// Equality of two spheres.
  template <class SphereT1, class RealT1, int DimN1, class SphereT2, class RealT2, int DimN2>
  inline bool operator==( SphereBase<SphereT1,RealT1,DimN1> const& sphere1, SphereBase<SphereT2,RealT2,DimN2> const& sphere2 ) {
    return sphere1.center_()==sphere2.center_() && sphere1.radius()==sphere2.radius();
  }

  /// Inequality of two spheres.
  template <class SphereT1, class RealT1, int DimN1, class SphereT2, class RealT2, int DimN2>
  inline bool operator!=( SphereBase<SphereT1,RealT1,DimN1> const& sphere1, SphereBase<SphereT2,RealT2,DimN2> const& sphere2 ) {
    return sphere1.center_()!=sphere2.center_() || sphere1.radius()!=sphere2.radius();
  }

  /// Writes a sphere to an ostream.
  template <class SphereT, class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, SphereBase<SphereT,RealT,DimN> const& sphere ) {
    return os << "(" << sphere.center_() << "-" << sphere.radius() << ")";
  }

  // *******************************************************************
  // class Sphere
  // *******************************************************************

  /// A general fixed-dimensional sphere class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT, int DimN = 0>
  class Sphere : public SphereBase<Sphere<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty sphere.
    Sphere() : SphereBase<Sphere<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a sphere with the given center and radius.
    template <class VectorT>
    Sphere( VectorBase<VectorT> const& center, RealT radius ) :
      SphereBase<Sphere<RealT, DimN>, RealT, DimN>( center, radius ) {}

    /// Constructs a 2D sphere with the given center point
    /// coordinates and radius.  (Only valid for 2D spheres.)
    Sphere( RealT centerx, RealT centery, RealT radius )
      : SphereBase<Sphere<RealT, DimN>, RealT, DimN>( centerx, centery, radius ) {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D sphere with the given center point
    /// coordinates and radius.  (Only valid for 3D spheres.)
    Sphere( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : SphereBase<Sphere<RealT, DimN>, RealT, DimN>( centerx, centery, centerz, radius ) {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Copy constructor.
    template <class SphereT1, class RealT1, int DimN1>
    Sphere( SphereBase<SphereT1, RealT1, DimN1> const& sphere )
      : SphereBase<Sphere<RealT, DimN>, RealT, DimN>( sphere.center_() , sphere.radius() ) {}

    /// Copy assignment operator.
    template <class SphereT1, class RealT1, int DimN1>
    Sphere& operator=( SphereBase<SphereT1, RealT1, DimN1> const& sphere ) {
      this->center_() = sphere.center_();
      this->radius() = sphere.radius();
      this->m_empty = sphere.empty();
      return *this;
    }
  };

  /// A general arbitrary-dimensional sphere class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT>
  class Sphere<RealT, 0> : public SphereBase<Sphere<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor. Constructs an empty sphere.
    Sphere() : SphereBase<Sphere<RealT, 0>, RealT, 0>() {}

    /// Constructs a sphere with the given center and radius.
    template <class VectorT>
    Sphere( VectorBase<VectorT> const& center, RealT radius ) :
      SphereBase<Sphere<RealT, 0>, RealT, 0>( center, radius ) {}

    /// Constructs a 2D sphere with the given center point
    /// coordinates and radius.  (Only valid for 2D spheres.)
    Sphere( RealT centerx, RealT centery, RealT radius )
      : SphereBase<Sphere<RealT, 0>, RealT, 0>( centerx, centery, radius ) {}

    /// Constructs a 3D sphere with the given center point
    /// coordinates and radius.  (Only valid for 3D spheres.)
    Sphere( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : SphereBase<Sphere<RealT, 0>, RealT, 0>( centerx, centery, centerz, radius ) {}

    /// Copy constructor.
    template <class SphereT1, class RealT1, int DimN1>
    Sphere( SphereBase<SphereT1, RealT1, DimN1> const& sphere )
      : SphereBase<Sphere<RealT, 0>, RealT, 0>( sphere.center_() , sphere.radius() ) {}

    /// Copy assignment operator.
    template <class SphereT1, class RealT1, int DimN1>
    Sphere& operator=( SphereBase<SphereT1, RealT1, DimN1> const& sphere ) {
      this->center_() = sphere.center_();
      this->radius() = sphere.radius();
      this->m_empty = sphere.empty();
      return *this;
    }
  };

} // namespace math

  // Convenience typedefs
  using geometry::Sphere;
  typedef Sphere<float64, 2> Sphere2;
  typedef Sphere<float64, 3> Sphere3;
  typedef Sphere<float64, 4> Sphere4;
  typedef Sphere<float64> SphereN;
  typedef Sphere<float32, 2> Sphere2f;
  typedef Sphere<float32, 3> Sphere3f;
  typedef Sphere<float32, 4> Sphere4f;
  typedef Sphere<float32> SphereNf;

} // namespace vw

#endif // __VW_GEOMETRY_SPHERE_H__
