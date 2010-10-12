// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BBox.h
///
/// Provides a generic bounding box.
#ifndef __VW_MATH_BBOX_H__
#define __VW_MATH_BBOX_H__

#include <iostream>
#include <limits>
#include <vector>
#include <cmath>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  /// \cond INTERNAL
  namespace vector_containment_comparison {
    template <class VectorT1, class VectorT2>
    inline bool operator<( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] < v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator<=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] <= v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator>( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] > v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator>=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] >= v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline VectorT1 max( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute max of vectors of different length." );
      VectorT1 v3;
      v3.set_size(v1.impl().size());
      for( size_t i=0; i<v1.impl().size(); ++i)
        v3[i] = std::max(v1.impl()[i], v2.impl()[i]);
      return v3;
    }

    template <class VectorT1, class VectorT2>
    inline VectorT1 min( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute min of vectors of different length." );
      VectorT1 v3;
      v3.set_size(v1.impl().size());
      for( size_t i=0; i<v1.impl().size(); ++i)
        v3[i] = std::min(v1.impl()[i], v2.impl()[i]);
      return v3;
    }

  } // namespace vector_containment_comparison
  /// \endcond


  // *******************************************************************
  // class BBoxBase
  // *******************************************************************

  /// A general n-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class BBoxT, class RealT, size_t DimN>
  class BBoxBase {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    BBoxBase() {
      // Make sure we have a type for which we know limits
      BOOST_STATIC_ASSERT(std::numeric_limits<RealT>::is_specialized);
      if (std::numeric_limits<RealT>::is_integer) {
        for (size_t i = 0; i < m_min.size(); ++i) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = std::numeric_limits<RealT>::min();
        }
      }
      else {
        for (size_t i = 0; i < m_min.size(); ++i) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = static_cast<RealT>(-std::numeric_limits<RealT>::max());
        }
      }
    }

    /// Constructs a bounding box with the given minimal and maximal
    /// points.
    template <class VectorT1, class VectorT2>
    BBoxBase( VectorBase<VectorT1> const& min, VectorBase<VectorT2> const& max ) :
      m_min( min ), m_max( max ) {}

    /// Returns the derived implementation type.
    BBoxT& impl() { return *static_cast<BBoxT*>(this); }

    /// Returns the derived implementation type.
    BBoxT const& impl() const { return *static_cast<BBoxT const*>(this); }

    /// Grows a bounding box to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
      if (m_min.size() == 0) {
        m_min = point;
        m_max = point;
      }
      else {
        for (size_t i = 0; i < m_min.size(); ++i) {
          if (point.impl()[i] > m_max[i])
            m_max[i] = RealT(point.impl()[i]);
          if (point.impl()[i] < m_min[i])
            m_min[i] = RealT(point.impl()[i]);
        }
      }
    }

    /// Grows a bounding box to include the given bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    void grow( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      grow(bbox.min());
      grow(bbox.max());
    }

    /// Crops (intersects) this bounding box to the given bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    void crop( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      for( size_t i=0; i<m_min.size(); ++i ) {
        if( m_min[i] < bbox.min()[i] ) {
          if( m_max[i] < bbox.min()[i] )
            m_min[i] = m_max[i];
          else
            m_min[i] = bbox.min()[i];
        }

        if( m_max[i] > bbox.max()[i] ) {
          if ( m_min[i] > bbox.max()[i] )
            m_max[i] = m_min[i];
          else
            m_max[i] = bbox.max()[i];
        }
      }
    }

    /// Expands this bounding box by the given offset in every direction.
    void expand( RealT offset ) {
      for( size_t i=0; i<m_min.size(); ++i ) {
        m_min[i] -= offset;
        m_max[i] += offset;
      }
    }

    /// Contracts this bounding box by the given offset in every direction.
    void contract( RealT offset ) {
      for( size_t i=0; i<m_min.size(); ++i ) {
        m_min[i] += offset;
        m_max[i] -= offset;
      }
    }

    /// Returns true if the given point is contained in the bounding box.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      using namespace vector_containment_comparison;
      VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
      return ((m_min.size() != 0) && (point >= m_min) && (point < m_max));
    }

    /// Returns true if the given bounding box is entirely contained
    /// in this bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    bool contains( const BBoxBase<BBoxT1, RealT1, DimN1> &bbox ) const {
      using namespace vector_containment_comparison;
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      return ((m_min.size() != 0) && (bbox.min() >= m_min) && (bbox.max() <= m_max));
    }

    /// Returns true if the given bounding box intersects this
    /// bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    bool intersects( const BBoxBase<BBoxT1, RealT1, DimN1>& bbox ) const {
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      for( size_t i=0; i<m_min.size(); ++i ) {
        if( m_min[i] >= bbox.max()[i] ||
            m_max[i] <= bbox.min()[i] )
          return false;
      }
      return (m_min.size() != 0);
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
      return impl().width();
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      return impl().height();
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      return impl().depth();
    }

    /// Returns true if the bounding box is empty (i.e. degenerate).
    bool empty() const {
      for( size_t i=0; i<m_min.size(); ++i )
        if( m_min[i] >= m_max[i] ) return true;
      return (m_min.size() <= 0);
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBoxBase& operator*=( ScalarT s ) {
      m_min *= s;
      m_max *= s;
      return *this;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBoxBase& operator/=( ScalarT s ) {
      m_min /= s;
      m_max /= s;
      return *this;
    }

    /// Offsets the bounding box by the given vector.
    template <class VectorT>
    BBoxBase& operator+=( VectorBase<VectorT> const& v ) {
      m_min += v;
      m_max += v;
      return *this;
    }

    /// Offsets the bounding box by the negation of the given vector.
    template <class VectorT>
    BBoxBase& operator-=( VectorBase<VectorT> const& v ) {
      m_min -= v;
      m_max -= v;
      return *this;
    }

  protected:
    Vector<RealT, DimN> m_min, m_max;
  };

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator*( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s ) {
    BBoxT result = bbox.impl();
    result *= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator/( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s ) {
    BBoxT result = bbox.impl();
    result /= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator*( ScalarT s, BBoxBase<BBoxT, RealT, DimN> const& bbox ) {
    return bbox * s;
  }

  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator+( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBoxT result = bbox.impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator+( VectorBase<VectorT> const& v, BBoxBase<BBoxT, RealT, DimN> const& bbox ) {
    return bbox + v;
  }

  /// Offsets a bounding box by the negation of the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator-( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBoxT result = bbox.impl();
    result -= v.impl();
    return result;
  }

  /// Equality of two bounding boxes.
  template <class BBoxT1, class RealT1, size_t DimN1, class BBoxT2, class RealT2, size_t DimN2>
  inline bool operator==( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 ) {
    return bbox1.min()==bbox2.min() && bbox1.max()==bbox2.max();
  }

  /// Inequality of two bounding boxes.
  template <class BBoxT1, class RealT1, size_t DimN1, class BBoxT2, class RealT2, size_t DimN2>
  inline bool operator!=( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 ) {
    return bbox1.min()!=bbox2.min() || bbox1.max()!=bbox2.max();
  }

  /// Writes a bounding box to an ostream.
  template <class BBoxT, class RealT, size_t DimN>
  std::ostream& operator<<( std::ostream& os, BBoxBase<BBoxT,RealT,DimN> const& bbox ) {
    return os << "(" << bbox.min() << "-" << bbox.max() << ")";
  }

  /// Asymmetricaly scale a bounding box, by elementwise vector product
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT elem_prod( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    return BBoxT( elem_prod(bbox.min(),v), elem_prod(bbox.max(),v) );
  }

  /// Asymmetricaly scale a bounding box, by elementwise vector quotient
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT elem_quot( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    return BBoxT( elem_quot(bbox.min(),v), elem_quot(bbox.max(),v) );
  }

  // *******************************************************************
  // class BBox
  // *******************************************************************

  /// A general fixed-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT, size_t DimN = 0>
  class BBox : public BBoxBase<BBox<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    BBox() : BBoxBase<BBox<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a bounding box with the given minimal and maximal
    /// points.
    BBox( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max ) :
      BBoxBase<BBox<RealT, DimN>, RealT, DimN>( min, max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>(
          Vector<RealT,2>(minx,miny),
          Vector<RealT,2>(minx+width,miny+height) ) {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>(
          Vector<RealT,3>(minx,miny,minz),
          Vector<RealT,3>(minx+width,miny+height,minz+depth) ) {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Copy constructor.
    template <class BBoxT1, class RealT1, size_t DimN1>
    BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>( bbox.min(), bbox.max() ) {}

    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, size_t DimN1>
    BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      this->min() = bbox.min();
      this->max() = bbox.max();
      return *this;
    }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      BOOST_STATIC_ASSERT( DimN >= 1 );
      return this->max()[0] - this->min()[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      BOOST_STATIC_ASSERT( DimN >= 2 );
      return this->max()[1] - this->min()[1];
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      BOOST_STATIC_ASSERT( DimN >= 3 );
      return this->max()[2] - this->min()[2];
    }
  };

  /// A general arbitrary-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT>
  class BBox<RealT, 0> : public BBoxBase<BBox<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    BBox() : BBoxBase<BBox<RealT, 0>, RealT, 0>() {}

    /// Constructs a bounding box with the given minimal and maximal
    /// points.
    BBox( Vector<RealT, 0> const& min, Vector<RealT, 0> const& max ) :
      BBoxBase<BBox<RealT, 0>, RealT, 0>( min, max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>(
          Vector<RealT,2>(minx,miny),
          Vector<RealT,2>(minx+width,miny+height) ) {}

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>(
          Vector<RealT,3>(minx,miny,minz),
          Vector<RealT,3>(minx+width,miny+height,minz+depth) ) {}

    /// Copy constructor.
    template <class BBoxT1, class RealT1, size_t DimN1>
    BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>( bbox.min(), bbox.max() ) {}

    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, size_t DimN1>
    BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      this->min() = bbox.min();
      this->max() = bbox.max();
      return *this;
    }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      VW_ASSERT(this->min().size() >= 1, LogicErr() << "BBox must be of dimension >= 1 to get width.");
      return this->max()[0] - this->min()[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      VW_ASSERT(this->min().size() >= 2, LogicErr() << "BBox must be of dimension >= 2 to get height.");
      return this->max()[1] - this->min()[1];
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      VW_ASSERT(this->min().size() >= 3, LogicErr() << "BBox must be of dimension >= 3 to get depth.");
      return this->max()[2] - this->min()[2];
    }

  };

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
  typedef BBox<float64> BBoxN;
  typedef BBox<float32> BBoxNf;
  typedef BBox<int32> BBoxNi;

} // namespace vw

#endif // __VW_MATH_BBOX_H__
