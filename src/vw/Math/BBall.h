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

/// \file BBall.h
///
/// Provides a generic bounding ball.
#ifndef __VW_MATH__BALL_H__
#define __VW_MATH__BALL_H__

#include <iostream>
#include <math.h>

#include <boost/static_assert.hpp>

#include <vw/Math/BShape.h>
#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  // *******************************************************************
  // class BBallBase
  // *******************************************************************

  /// A general n-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class BBallT, class RealT, int DimN>
  class BBallBase : public BShapeBase<BBallBase<BBallT, RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBallBase() : m_radius( 0 ) {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBallBase( VectorBase<VectorT> const& center, RealT radius ) :
      m_center( center ), m_radius( radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBallBase( RealT centerx, RealT centery, RealT radius )
      : m_center( Vector<RealT,2>(centerx,centery) ), 
        m_radius( radius ) {}

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBallBase( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : m_center( Vector<RealT,3>(centerx,centery,centerz) ), 
        m_radius( radius ) {}

    /// Returns the derived implementation type.
    BBallT& impl() { return *static_cast<BBallT*>(this); }
    
    /// Returns the derived implementation type.
    BBallT const& impl() const { return *static_cast<BBallT const*>(this); }

    /// Grows a bounding ball to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      Vector<RealT, DimN> v = point - m_center;
      RealT dist = norm_2(v);
      RealT radius_d;
      if (dist > m_radius) {
        radius_d = (dist - m_radius)/2;
        m_center += v/dist*radius_d;
        m_radius += radius_d;
      }
    }
    
    /// Grows a bounding ball to include the given bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    void grow( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      Vector<RealT, DimN> v = bball.m_center - m_center;
      grow(bball.m_center + normalize(v)*bball.m_radius);
    }

    /// Crops (intersects) this bounding ball to the given bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    void crop( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      Vector<RealT, DimN> v = bball.m_center - m_center;
      RealT dist = norm_2(v);
      if (dist > (m_radius + bball.m_radius)) {
        // intersects(bball) == false
        m_radius = 0;
      }
      else if (dist <= (m_radius - bball.m_radius)) {
        // contains(bball) == true
        m_center = bball.m_center;
        m_radius = bball.m_radius;
      }
      else if (dist <= (bball.m_radius - m_radius)) {
        // bball.contains(*this) == true
        // smallest BBall is *this
        return;
      }
      else if (m_radius >= dist) {
        // smallest BBall is bball
        m_center = bball.m_center;
        m_radius = bball.m_radius;
      }
      else if (bball.m_radius >= dist) {
        // smallest BBall is *this
        return;
      }
      else {
        // normal intersection
        RealT x = (dist + (m_radius*m_radius - bball.m_radius*bball.m_radius)/dist)/2;
        m_center += v/dist*x;
        m_radius = std::sqrt(m_radius*m_radius - x*x);
      }
    }

    /// Expands this bounding ball by the given offset in every direction.
    void expand( RealT offset ) {
      m_radius += offset;
    }

    /// Contracts this bounding ball by the given offset in every direction.
    void contract( RealT offset ) {
      m_radius -= offset;
    }

    /// Returns true if the given point is contained in the bounding ball.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      VW_ASSERT(point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      return (norm_2(m_center - point) <= m_radius);
    }

    /// Returns true if the given bounding ball is entirely contained
    /// in this bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    bool contains( const BBallBase<BBallT1, RealT1, DimN1> &bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius - bball.m_radius));
    }

    /// Returns true if the given bounding ball intersects this
    /// bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    bool intersects( const BBallBase<BBallT1, RealT1, DimN1>& bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius + bball.m_radius));
    }

    /// Returns the size (i.e. the diameter) of the bounding ball.
    RealT size() const { return 2*m_radius; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> center() const { return m_center; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> const& center_() const { return m_center; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> &center_() { return m_center; }

    /// Returns the radius of the bounding ball.
    Vector<RealT, DimN> const& radius() const { return m_radius; }

    /// Returns the radius of the bounding ball.
    Vector<RealT, DimN> &radius() { return m_radius; }

    /// Returns true if the bounding ball is empty (i.e. degenerate).
    bool empty() const {
      return (m_radius <= 0);
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBallBase& operator*=( ScalarT s ) {
      m_center *= s;
      m_radius *= s;
      return *this;
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBallBase& operator/=( ScalarT s ) {
      m_center /= s;
      m_radius /= s;
      return *this;
    }

    /// Offsets the bounding ball by the given vector.
    template <class VectorT>
    BBallBase& operator+=( VectorBase<VectorT> const& v ) {
      m_center += v;
      return *this;
    }

    /// Offsets the bounding ball by the negation of the given vector.
    template <class VectorT>
    BBallBase& operator-=( VectorBase<VectorT> const& v ) {
      m_center -= v;
      return *this;
    }

  private:
    Vector<RealT, DimN> m_center;
    RealT m_radius;
  };
  
  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator*( BBallBase<BBallT, RealT, DimN> const& bball, ScalarT s ) {
    BBallT result = bball.impl();
    result *= s;
    return result;
  }

  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator/( BBallBase<BBallT, RealT, DimN> const& bball, ScalarT s ) {
    BBallT result = bball.impl();
    result /= s;
    return result;
  }

  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator*( ScalarT s, BBallBase<BBallT, RealT, DimN> const& bball ) {
    return bball * s;
  }
  
  /// Offsets a bounding ball by the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator+( BBallBase<BBallT, RealT, DimN> const& bball, VectorBase<VectorT> const& v ) {
    BBallT result = bball.impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding ball by the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator+( VectorBase<VectorT> const& v, BBallBase<BBallT, RealT, DimN> const& bball ) {
    return bball + v;
  }

  /// Offsets a bounding ball by the negation of the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator-( BBallBase<BBallT, RealT, DimN> const& bball, VectorBase<VectorT> const& v ) {
    BBallT result = bball.impl();
    result -= v.impl();
    return result;
  }
  
  /// Equality of two bounding balls.
  template <class BBallT1, class RealT1, int DimN1, class BBallT2, class RealT2, int DimN2>
  inline bool operator==( BBallBase<BBallT1,RealT1,DimN1> const& bball1, BBallBase<BBallT2,RealT2,DimN2> const& bball2 ) {
    return bball1.center_()==bball2.center_() && bball1.radius()==bball2.radius();
  }
  
  /// Inequality of two bounding balls.
  template <class BBallT1, class RealT1, int DimN1, class BBallT2, class RealT2, int DimN2>
  inline bool operator!=( BBallBase<BBallT1,RealT1,DimN1> const& bball1, BBallBase<BBallT2,RealT2,DimN2> const& bball2 ) {
    return bball1.center_()!=bball2.center_() || bball1.radius()!=bball2.radius();
  }
  
  /// Writes a bounding ball to an ostream.
  template <class BBallT, class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBallBase<BBallT,RealT,DimN> const& bball ) {
    return os << "(" << bball.center_() << "-" << bball.radius() << ")";
  }

  // *******************************************************************
  // class BBall
  // *******************************************************************

  /// A general fixed-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT, int DimN = 0>
  class BBall : public BBallBase<BBall<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBall() : BBallBase<BBall<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBall( VectorBase<VectorT> const& center, RealT radius ) :
      BBallBase<BBall<RealT, DimN>, RealT, DimN>( center, radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT radius )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( centerx, centery, radius ) {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( centerx, centery, centerz, radius ) {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Copy constructor.
    template <class BBallT1, class RealT1, int DimN1>
    BBall( BBallBase<BBallT1, RealT1, DimN1> const& bball )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( bball.center_() , bball.radius() ) {}
    
    /// Copy assignment operator.
    template <class BBallT1, class RealT1, int DimN1>
    BBall& operator=( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      this->center_() = bball.center_();
      this->radius() = bball.radius();
      return *this;
    }
  };

  /// A general arbitrary-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT>
  class BBall<RealT, 0> : public BBallBase<BBall<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBall() : BBallBase<BBall<RealT, 0>, RealT, 0>() {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBall( VectorBase<VectorT> const& center, RealT radius ) :
      BBallBase<BBall<RealT, 0>, RealT, 0>( center, radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT radius )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( centerx, centery, radius ) {}

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( centerx, centery, centerz, radius ) {}

    /// Copy constructor.
    template <class BBallT1, class RealT1, int DimN1>
    BBall( BBallBase<BBallT1, RealT1, DimN1> const& bball )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( bball.center_() , bball.radius() ) {}
    
    /// Copy assignment operator.
    template <class BBallT1, class RealT1, int DimN1>
    BBall& operator=( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      this->center_() = bball.center_();
      this->radius() = bball.radius();
      return *this;
    }
  };

} // namespace math

  // Convenience typedefs
  using math::BBall;
  typedef BBall<float64, 2> BBall2;
  typedef BBall<float64, 3> BBall3;
  typedef BBall<float64, 4> BBall4;
  typedef BBall<float32, 2> BBall2f;
  typedef BBall<float32, 3> BBall3f;
  typedef BBall<float32, 4> BBall4f;
  typedef BBall<int32, 2> BBall2i;
  typedef BBall<int32, 3> BBall3i;
  typedef BBall<int32, 4> BBall4i;
  typedef BBall<float64> BBallN;
  typedef BBall<float32> BBallNf;
  typedef BBall<int32> BBallNi;

} // namespace vw

#endif // __VW_MATH__BALL_H__
