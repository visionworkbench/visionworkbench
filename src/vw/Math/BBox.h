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


/// \file BBox.h
///
/// Provides a generic bounding box.
#ifndef __VW_MATH_BBOX_H__
#define __VW_MATH_BBOX_H__

#include <vw/Math/Vector.h>

#include <iostream>
#include <limits>
#include <vector>
#include <cmath>

#include <boost/static_assert.hpp>

namespace vw {
namespace math {

  // TODO: Move these to Vector.h and delete this namespace
  /// \cond INTERNAL
  namespace vector_containment_comparison {

    /// Implement < operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline bool operator<( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different lengths." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] < v2.impl()[i] ) ) return false;
      return true;
    }

    /// Implement <= operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline bool operator<=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different lengths." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] <= v2.impl()[i] ) ) return false;
      return true;
    }

    /// Implement > operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline bool operator>( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different lengths." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] > v2.impl()[i] ) ) return false;
      return true;
    }

    /// Implement >= operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline bool operator>=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different lengths." );
      for( size_t i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] >= v2.impl()[i] ) ) return false;
      return true;
    }

    /// Implement max() operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline VectorT1 max( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute max of vectors of different lengths." );
      VectorT1 v3;
      v3.set_size(v1.impl().size());
      for( size_t i=0; i<v1.impl().size(); ++i)
        v3[i] = std::max(v1.impl()[i], v2.impl()[i]);
      return v3;
    }

    /// Implement min() operation for VectorBase
    template <class VectorT1, class VectorT2>
    inline VectorT1 min( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute min of vectors of different lengths." );
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
  /// represented by vectors pointing to the minimal and maximal corners.
  template <class BBoxT, class RealT, size_t DimN>
  class BBoxBase {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    inline BBoxBase();

    /// Constructs a bounding box with the given minimal and maximal points.
    template <class VectorT1, class VectorT2>
    BBoxBase( VectorBase<VectorT1> const& min, VectorBase<VectorT2> const& max ) : m_min( min ), m_max( max ) {}

    BBoxT      & impl()       { return *static_cast<BBoxT      *>(this); } ///< Returns the derived implementation type.
    BBoxT const& impl() const { return *static_cast<BBoxT const*>(this); } ///< Returns the derived implementation type.

    /// Returns true if the bounding box is empty (i.e. degenerate).
    inline bool empty() const;
    
    inline Vector<RealT, DimN>        size  () const; ///< Returns the size (i.e. the diagonal vector) of the bounding box.
    inline Vector<RealT, DimN>        center() const;                  ///< Returns the center point of the bounding box.
    Vector<RealT, DimN> const& min   () const { return m_min; } ///< Returns the minimal point of the bounding box.
    Vector<RealT, DimN>      & min   ()       { return m_min; } ///< Returns the minimal point of the bounding box.
    Vector<RealT, DimN> const& max   () const { return m_max; } ///< Returns the maximal point of the bounding box.   
    Vector<RealT, DimN>      & max   ()       { return m_max; } ///< Returns the maximal point of the bounding box.

    inline RealT width () const; ///< Returns the width  (i.e. size in the first  dimension) of the bounding box.
    inline RealT height() const; ///< Returns the height (i.e. size in the second dimension) of the bounding box.
    inline RealT depth () const; ///< Returns the depth  (i.e. size in the third  dimension) of the bounding box.

    inline void expand  ( RealT offset ); ///< Expands   this bounding box by the given offset in every direction.
    inline void contract( RealT offset ); ///< Contracts this bounding box by the given offset in every direction.


    /// Grows a bounding box to include the given point.
    template <class VectorT>
    inline void grow( VectorBase<VectorT> const& point );

    /// Grows a bounding box to include the given bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline void grow( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) { grow(bbox.min()); grow(bbox.max()); }

    /// Crops (intersects) this bounding box to the given bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline void crop( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

    /// Returns true if the given point is contained in the bounding box.
    template <class VectorT>
    inline bool contains( const VectorBase<VectorT> &point ) const;

    /// Returns true if the given bounding box is entirely contained in this bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline bool contains( const BBoxBase<BBoxT1, RealT1, DimN1> &bbox ) const;

    /// Returns true if the given bounding box intersects this bounding box.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline bool intersects( const BBoxBase<BBoxT1, RealT1, DimN1>& bbox ) const;


    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    inline BBoxBase& operator*=( ScalarT s );

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    inline BBoxBase& operator/=( ScalarT s );

    /// Offsets the bounding box by the given vector.
    template <class VectorT>
    inline BBoxBase& operator+=( VectorBase<VectorT> const& v );

    /// Offsets the bounding box by the negation of the given vector.
    template <class VectorT>
    inline BBoxBase& operator-=( VectorBase<VectorT> const& v );

  protected:
    Vector<RealT, DimN> m_min, m_max;
  }; // End class BBoxBase

  // Overloaded functions for operating on BBoxBase
  // - Some of these operations should probably have been named class functions.

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator*( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s );

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator/( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s );

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, size_t DimN, class ScalarT>
  inline BBoxT operator*( ScalarT s, BBoxBase<BBoxT, RealT, DimN> const& bbox ) { return bbox * s; }

  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator+( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v );

  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator+( VectorBase<VectorT> const& v, BBoxBase<BBoxT, RealT, DimN> const& bbox ) { return bbox + v; }

  /// Offsets a bounding box by the negation of the given vector.
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT operator-( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v );

  /// Equality of two bounding boxes.
  template <class BBoxT1, class RealT1, size_t DimN1, class BBoxT2, class RealT2, size_t DimN2>
  inline bool operator==( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 );

  /// Inequality of two bounding boxes.
  template <class BBoxT1, class RealT1, size_t DimN1, class BBoxT2, class RealT2, size_t DimN2>
  inline bool operator!=( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 );

  /// Writes a bounding box to an ostream.
  template <class BBoxT, class RealT, size_t DimN>
  std::ostream& operator<<( std::ostream& os, BBoxBase<BBoxT,RealT,DimN> const& bbox ) {
    return os << "(" << bbox.min() << "-" << bbox.max() << ")";
  }
  

  /// Specialization for 2d boxes
  template <class BBoxT, class RealT>
  std::ostream& operator<<( std::ostream& os, BBoxBase<BBoxT,RealT,2> const& bbox );

  /// Asymmetricaly scale a bounding box, by elementwise vector product
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT elem_prod( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v );

  /// Asymmetricaly scale a bounding box, by elementwise vector quotient
  template <class BBoxT, class RealT, size_t DimN, class VectorT>
  inline BBoxT elem_quot( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v );

  // End functions for BBoxBase
  //--------------------------------------------------------------------

  // *******************************************************************
  // class BBox - Fixed dimensions
  // *******************************************************************

  /// A general fixed-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal corners.
  template <class RealT, size_t DimN = 0>
  class BBox : public BBoxBase<BBox<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    inline BBox() : BBoxBase<BBox<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a bounding box with the given minimal and maximal points.
    inline BBox( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max ) :
      BBoxBase<BBox<RealT, DimN>, RealT, DimN>( min, max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bounding boxes.)
    inline BBox( RealT minx, RealT miny, RealT width, RealT height );

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bouding boxes.)
    inline BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth );

    /// Copy constructor.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

    /// Returns true if the bounding box is empty (i.e. degenerate).
    inline bool empty() const;
    
    inline RealT width () const; ///< Returns the width  (i.e. size in the first  dimension) of the bounding box.   
    inline RealT height() const; ///< Returns the height (i.e. size in the second dimension) of the bounding box.
    inline RealT depth () const; ///< Returns the depth  (i.e. size in the third  dimension) of the bounding box.
  };

  // *******************************************************************
  // class BBox - Variable dimensions
  // *******************************************************************

  /// A general arbitrary-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal corners.
  template <class RealT>
  class BBox<RealT, 0> : public BBoxBase<BBox<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding
    /// box, whose limits are at the opposite corners of the underlying
    /// numeric space.  This is a useful starting point if you intend
    /// to grow your bounding box to fit a collection of items.
    inline BBox() : BBoxBase<BBox<RealT, 0>, RealT, 0>() {}

    /// Constructs a bounding box with the given minimal and maximal points.
    inline BBox( Vector<RealT, 0> const& min, Vector<RealT, 0> const& max );

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bounding boxes.)
    inline BBox( RealT minx, RealT miny, RealT width, RealT height );

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bounding boxes.)
    inline BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth );

    /// Copy constructor.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, size_t DimN1>
    inline BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

    /// Returns true if the bounding box is empty (i.e. degenerate).
    inline bool empty() const;

    inline RealT width () const; ///< Returns the width  (i.e. size in the first  dimension) of the bounding box.   
    inline RealT height() const; ///< Returns the height (i.e. size in the second dimension) of the bounding box.
    inline RealT depth () const; ///< Returns the depth  (i.e. size in the third  dimension) of the bounding box.
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
  typedef BBox<int32,   2> BBox2i;
  typedef BBox<int32,   3> BBox3i;
  typedef BBox<int32,   4> BBox4i;
  typedef BBox<uint32,  2> BBox2u;
  typedef BBox<uint32,  3> BBox3u;
  typedef BBox<uint32,  4> BBox4u;
  typedef BBox<float64   > BBoxN;
  typedef BBox<float32   > BBoxNf;
  typedef BBox<int32     > BBoxNi;

  /// A helper function to grow a floating-point bounding box
  /// to the smallest enclosing integer bounding box.
  template <class BBoxT1, class RealT1, size_t DimN1>
  inline typename boost::enable_if<boost::is_float<RealT1>,BBox<int32,DimN1> >::type
  grow_bbox_to_int( math::BBoxBase<BBoxT1, RealT1, DimN1> const& bbox );

  // For BBoxes that are already int, do nothing
  template <class BBoxT1, class RealT1, size_t DimN1>
  inline typename boost::enable_if<boost::is_integral<RealT1>,BBox<RealT1,DimN1> >::type
  grow_bbox_to_int( math::BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) { return bbox; }
} // namespace vw


#include "BBox.tcc"

#endif // __VW_MATH_BBOX_H__
