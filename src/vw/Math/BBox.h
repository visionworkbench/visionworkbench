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

/// \file BBox.h
///
/// Provides generic bounding shapes.
#ifndef __VW_MATH__BBOX_H__
#define __VW_MATH__BBOX_H__

#include <limits>
#include <math.h>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  /// \cond INTERNAL
  namespace vector_containment_comparison {
    template <class ElemT, int SizeN>
    inline bool operator<( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] < v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator<=( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] <= v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator>( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] > v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator>=( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] >= v2[i] ) ) return false;
      return true;
    }

  } // namespace vector_containment_comparison
  /// \endcond


  // *******************************************************************
  // class BShapeBase
  // *******************************************************************

  /// A CRTP base class for general n-dimensional bounding shapes.  
  /// Provides a mechanism for restricting function arguments to 
  /// bounding shapes, provides general bounding-shape operations,
  /// and provides the various arithmetic assignment operators.
  template <class BShapeT, class BShapeRealT, int BShapeDimN>
  class BShapeBase {
  typedef Vector<BShapeRealT, BShapeDimN> BShapeVectorT;
  public:

    /// Returns the derived implementation type.
    BShapeT& shape_impl() { return *static_cast<BShapeT*>(this); }
    
    /// Returns the derived implementation type.
    BShapeT const& shape_impl() const { return *static_cast<BShapeT const*>(this); }

    /// Grows a bounding shape to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      shape_impl().grow(point);
    }
    
    /// Grows a bounding shape to include the given bounding shape.
    void grow( BShapeT const& bshape ) {
      shape_impl().grow(bshape);
    }
    
    /// Crops (intersects) this bounding shape to the given bounding shape.
    void crop( BShapeT const& bshape ) {
      shape_impl().crop(bshape);
    }

    /// Expands this bounding shape by the given offset in every direction.
    void expand( BShapeRealT offset ) {
      shape_impl().expand(offset);
    }

    /// Contracts this bounding shape by the given offset in every direction.
    void contract( BShapeRealT offset ) {
      shape_impl().contract(offset);
    }

    /// Returns true if the given point is contained in the bounding shape.
    template <class VectorT>
    bool contains( VectorBase<VectorT> const& point ) const {
      return shape_impl().contains(point);
    }
    
    /// Returns true if the given bounding shape is entirely contained
    /// in this bounding shape.
    bool contains( const BShapeT &bshape ) const {
      return shape_impl().contains(bshape);
    }
    
    /// Returns true if the given bounding shape intersects this
    /// bounding shape.
    bool intersects( const BShapeT& bshape ) const {
      return shape_impl().intersects(bshape);
    }

    /// Returns the center point of the bounding shape.
    BShapeVectorT center() const {
      return shape_impl().center();
    }

    /// Returns true if the bounding shape is empty (i.e. degenerate).
    bool empty() const {
      return shape_impl().empty();
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    BShapeT& operator*=( ScalarT s ) {
      return shape_impl() = shape_impl() * s;
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    BShapeT& operator/=( ScalarT s ) {
      return shape_impl() = shape_impl() / s;
    }

    /// Offsets the bounding shape by the given vector.
    template <class VectorT>
    BShapeT& operator+=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() + v;
    }

    /// Offsets the bounding shape by the negation of the given vector.
    template <class VectorT>
    BShapeT& operator-=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() - v;
    }
  };
  
  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class ScalarT>
  inline BShapeT operator*( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape, ScalarT s ) {
    BShapeT result = bshape.shape_impl();
    result *= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class ScalarT>
  inline BShapeT operator/( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape, ScalarT s ) {
    BShapeT result = bshape.shape_impl();
    result /= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class ScalarT>
  inline BShapeT operator*( ScalarT s, BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape ) {
    return bshape * s;
  }
  
  /// Offsets a bounding shape by the given vector.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class VectorT>
  inline BShapeT operator+( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape, VectorBase<VectorT> const& v ) {
    BShapeT result = bshape.shape_impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding shape by the given vector.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class VectorT>
  inline BShapeT operator+( VectorBase<VectorT> const& v, BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape ) {
    return bshape + v;
  }

  /// Offsets a bounding shape by the negation of the given vector.
  template <class BShapeT, class BShapeRealT, int BShapeDimN, class VectorT>
  inline BShapeT operator-( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape, VectorBase<VectorT> const& v ) {
    BShapeT result = bshape.shape_impl();
    result -= v.impl();
    return result;
  }
  
  /// Equality of two bounding shapes.
  template <class BShapeT, class BShapeRealT, int BShapeDimN>
  inline bool operator==( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape1, BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape2 ) {
    return bshape1.shape_impl() == bshape2.shape_impl();
  }
  
  /// Inequality of two bounding shapes.
  template <class BShapeT, class BShapeRealT, int BShapeDimN>
  inline bool operator!=( BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape1, BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape2 ) {
    return bshape1.shape_impl() != bshape2.shape_impl();
  }

  /// Writes a bounding shape to an ostream.
  template <class BShapeT, class BShapeRealT, int BShapeDimN>
  std::ostream& operator<<( std::ostream& os, BShapeBase<BShapeT, BShapeRealT, BShapeDimN> const& bshape ) {
    return os << bshape.shape_impl();
  }

  // *******************************************************************
  // class BBox
  // *******************************************************************

  /// A general n-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT, int DimN>
  class BBox : public BShapeBase<BBox<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding 
    /// box, whose limits are at the opposite corners of the underlying 
    /// numeric space.  This is a useful starting point if you intend 
    /// to grow your bounding box to fit a collection of items.
    BBox() {
      // Make sure we have a type for which we know limits
      BOOST_STATIC_ASSERT(std::numeric_limits<RealT>::is_specialized);
      if (std::numeric_limits<RealT>::is_integer) {
        for (int i = 0; i < DimN; i++) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = std::numeric_limits<RealT>::min();
        }
      }
      else {
        for (int i = 0; i < DimN; i++) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = -std::numeric_limits<RealT>::max();
        }
      }
    }

    /// Constructs a bounding box with the given minimal and maximal 
    /// points.
    BBox( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max ) :
      m_min( min ), m_max( max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : m_min( Vector<RealT,2>(minx,miny) ), 
        m_max( Vector<RealT,2>(minx+width,miny+height) )
    {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Standard copy constructor.
    BBox( BBox const& bbox ) : m_min( bbox.m_min ), m_max( bbox.m_max ) {}
    
    /// Standard copy assignment operator.
    BBox& operator=( BBox const& bbox ) {
      m_min = bbox.m_min;
      m_max = bbox.m_max;
      return *this;
    }

    /// Grows a bounding box to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(point.impl().size() == DimN, ArgumentErr() << "Vector must have dimension " << DimN << ".");
      for (int i = 0; i < DimN; i++) {
	if (point.impl()[i] > m_max[i])
	  m_max[i] = point.impl()[i];
	if (point.impl()[i] < m_min[i])
	  m_min[i] = point.impl()[i];
      }
    }
    
    /// Grows a bounding box to include the given bounding box.
    void grow( BBox const& bbox ) {
      grow(bbox.min());
      grow(bbox.max());
    }

    /// Crops (intersects) this bounding box to the given bounding box.
    void crop( BBox const& bbox ) {
      for( int i=0; i<DimN; ++i ) {
        if( m_min[i] < bbox.m_min[i] )
          if( m_max[i] < bbox.m_min[i] )
            m_min[i] = m_max[i]; 
          else
            m_min[i] = bbox.m_min[i];

        if( m_max[i] > bbox.m_max[i] )
          if ( m_min[i] > bbox.m_max[i] )
            m_max[i] = m_min[i];
          else
            m_max[i] = bbox.m_max[i];
      }
    }

    /// Expands this bounding box by the given offset in every direction.
    void expand( RealT offset ) {
      for( int i=0; i<DimN; ++i ) {
        m_min[i] -= offset;
        m_max[i] += offset;
      }
    }

    /// Contracts this bounding box by the given offset in every direction.
    void contract( RealT offset ) {
      for( int i=0; i<DimN; ++i ) {
        m_min[i] += offset;
        m_max[i] -= offset;
      }
    }

    /// Returns true if the given point is contained in the bounding box.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      using namespace vector_containment_comparison;
      VW_ASSERT(point.impl().size() == DimN, ArgumentErr() << "Vector must have dimension " << DimN << ".");
      Vector<RealT,DimN> point_ = point;
      return ((point_ >= m_min) && (point_ < m_max));
    }

    /// Returns true if the given bounding box is entirely contained
    /// in this bounding box.
    bool contains( const BBox &bbox ) const {
      using namespace vector_containment_comparison;
      return ((bbox.m_min >= m_min) && (bbox.m_max <= m_max));
    }

    /// Returns true if the given bounding box intersects this
    /// bounding box.
    bool intersects( const BBox& bbox ) const {
      for( int i=0; i<DimN; ++i ) {
        if( m_min[i] >= bbox.m_max[i] ||
            m_max[i] <= bbox.m_min[i] )
          return false;
      }
      return true;
    }

    /// Returns the size (i.e. the diagonal vector) of the bounding box.
    Vector<RealT, DimN> size() const { return (m_max - m_min); }

    /// Returns the center point of the bounding box.
    Vector<RealT, DimN> center() const { return 0.5 * (m_min + m_max); }

    /// Returns the minimal point of the bounding box.
    Vector<RealT, DimN> const& min() const { return m_min; }

    /// Returns the maximal point of the bounding box.
    Vector<RealT, DimN> const& max() const { return m_max; }

    /// Returns the minimal point of the bounding box (non-const overload).
    Vector<RealT, DimN> &min() { return m_min; }

    /// Returns the maximal point of the bounding box (non-const overload).
    Vector<RealT, DimN> &max() { return m_max; }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      BOOST_STATIC_ASSERT( DimN >= 1 );
      return m_max[0] - m_min[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      BOOST_STATIC_ASSERT( DimN >= 2 );
      return m_max[1] - m_min[1];
    }

    /// Returns true if the bounding box is empty (i.e. degenerate).
    bool empty() const {
      for( int i=0; i<DimN; ++i )
        if( m_min[i] >= m_max[i] ) return true;
      return false;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBox& operator*=( ScalarT s ) {
      m_min *= s;
      m_max *= s;
      return *this;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBox& operator/=( ScalarT s ) {
      m_min /= s;
      m_max /= s;
      return *this;
    }

    /// Offsets the bounding box by the given vector.
    template <class VectorT>
    BBox& operator+=( VectorBase<VectorT> const& v ) {
      m_min += v;
      m_max += v;
      return *this;
    }

    /// Offsets the bounding box by the negation of the given vector.
    template <class VectorT>
    BBox& operator-=( VectorBase<VectorT> const& v ) {
      m_min -= v;
      m_max -= v;
      return *this;
    }

  private:
    Vector<RealT, DimN> m_min, m_max;
  };

  /// Equality of two bounding boxes.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator==( BBox<Real1T,DimN> const& bbox1, BBox<Real2T,DimN> const& bbox2 ) {
    return bbox1.min()==bbox2.min() && bbox1.max()==bbox2.max();
  }

  /// Inequality of two bounding boxes.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator!=( BBox<Real1T,DimN> const& bbox1, BBox<Real2T,DimN> const& bbox2 ) {
    return bbox1.min()!=bbox2.min() || bbox1.max()!=bbox2.max();
  }
  
  /// Writes a bounding box to an ostream.
  template <class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBox<RealT,DimN> const& bbox ) {
    return os << "(" << bbox.min() << "-" << bbox.max() << ")";
  }

  // *******************************************************************
  // class BBall
  // *******************************************************************

  /// A general n-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT, int DimN>
  class BBall : public BShapeBase<BBall<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBall() : m_radius( 0 ) {}

    /// Constructs a bounding ball with the given center and radius.
    BBall( Vector<RealT, DimN> const& center, RealT radius ) :
      m_center( center ), m_radius( radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT radius )
      : m_center( Vector<RealT,2>(centerx,centery) ), 
        m_radius( radius )
    {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Standard copy constructor.
    BBall( BBall const& bball ) : m_center( bball.m_center ), m_radius( bball.m_radius ) {}
    
    /// Standard copy assignment operator.
    BBall& operator=( BBall const& bball ) {
      m_center = bball.m_center;
      m_radius = bball.m_radius;
      return *this;
    }

    /// Grows a bounding ball to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(point.impl().size() == DimN, ArgumentErr() << "Vector must have dimension " << DimN << ".");
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
    void grow( BBall const& bball ) {
      Vector<RealT, DimN> v = bball.m_center - m_center;
      grow(bball.m_center + normalize(v)*bball.m_radius);
    }

    /// Crops (intersects) this bounding ball to the given bounding ball.
    void crop( BBall const& bball ) {
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
      VW_ASSERT(point.impl().size() == DimN, ArgumentErr() << "Vector must have dimension " << DimN << ".");
      return (norm_2(m_center - point) <= m_radius);
    }

    /// Returns true if the given bounding ball is entirely contained
    /// in this bounding ball.
    bool contains( const BBall &bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius - bball.m_radius));
    }

    /// Returns true if the given bounding ball intersects this
    /// bounding ball.
    bool intersects( const BBall& bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius + bball.m_radius));
    }

    /// Returns the size (i.e. the diameter) of the bounding ball.
    RealT size() const { return 2*m_radius; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> center() const { return m_center; }

    /// Returns the radius of the bounding ball.
    Vector<RealT, DimN> radius() const { return m_radius; }

    /// Returns true if the bounding ball is empty (i.e. degenerate).
    bool empty() const {
      return (m_radius <= 0);
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBall& operator*=( ScalarT s ) {
      m_center *= s;
      m_radius *= s;
      return *this;
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBall& operator/=( ScalarT s ) {
      m_center /= s;
      m_radius /= s;
      return *this;
    }

    /// Offsets the bounding ball by the given vector.
    template <class VectorT>
    BBall& operator+=( VectorBase<VectorT> const& v ) {
      m_center += v;
      return *this;
    }

    /// Offsets the bounding ball by the negation of the given vector.
    template <class VectorT>
    BBall& operator-=( VectorBase<VectorT> const& v ) {
      m_center -= v;
      return *this;
    }

  private:
    Vector<RealT, DimN> m_center;
    RealT m_radius;
  };
  
  /// Equality of two bounding balls.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator==( BBall<Real1T,DimN> const& bball1, BBall<Real2T,DimN> const& bball2 ) {
    return bball1.center()==bball2.center() && bball1.radius()==bball2.radius();
  }
  
  /// Inequality of two bounding balls.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator!=( BBall<Real1T,DimN> const& bball1, BBall<Real2T,DimN> const& bball2 ) {
    return bball1.center()!=bball2.center() || bball1.radius()!=bball2.radius();
  }
  
  /// Writes a bounding ball to an ostream.
  template <class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBall<RealT,DimN> const& bball ) {
    return os << "(" << bball.center() << "-" << bball.radius() << ")";
  }

} // namespace math

  // Convenience typedefs
  using math::BBox;
  typedef BBox<float64, 2> BBox2;
  typedef BBox<float64, 3> BBox3;
  typedef BBox<float64, 4> BBox4;
  typedef BBox<float32, 2> BBox2f;
  typedef BBox<float32, 3> BBox3f;
  typedef BBox<float32, 4> BBox4f;
  typedef BBox<int32, 2> BBox2i;
  typedef BBox<int32, 3> BBox3i;
  typedef BBox<int32, 4> BBox4i;

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
} // namespace vw

#endif // __VW_MATH__BBOX_H__
