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

/// \file BConvex.h
///
/// Provides a generic convex bounding shape.
#ifndef __VW_MATH__BCONVEX_H__
#define __VW_MATH__BCONVEX_H__

#include <iostream>
#include <math.h>
#include <vector>

#include <boost/static_assert.hpp>

#include <vw/Math/BShape.h>
#include <vw/Math/Vector.h>

#include <ppl.hh> // Parma Polyhedra Library

namespace vw {
namespace math {
  
  // *******************************************************************
  // class BConvex
  // *******************************************************************
  
  /// A general arbitrary-dimensional convex shape class.
  template <class RealT>
  class BConvex : public BShapeBase<BConvex<RealT>, RealT, 0> {
  public:
    /// Dimension-only constructor.
    BConvex( unsigned dim ) : m_poly( dim, Parma_Polyhedra_Library::EMPTY ) {}
    
    /// Vector-of-points constructor.
    template <class ContainerT>
    BConvex( std::vector<ContainerT> const& points ) : m_poly( points[0].size(), Parma_Polyhedra_Library::EMPTY ) {
      unsigned num_points = points.size();
      for (unsigned i = 0; i < num_points; i++)
        grow(points[i]);
    }

    /// Destructor.
    ~BConvex() {}
    
    /// Copy constructor.
    template <class RealT1>
    BConvex( BConvex<RealT1> const& bconv )
      : m_poly( bconv.poly() ) {}
    
    /// Copy assignment operator.
    template <class RealT1>
    BConvex& operator=( BConvex<RealT1> const& bconv ) {
      m_poly = bconv.poly();
      return *this;
    }
    
    /// Returns the internal polyhedron.
    Parma_Polyhedra_Library::C_Polyhedron &poly() {return m_poly;}
    
    /// Returns the internal polyhedron.
    Parma_Polyhedra_Library::C_Polyhedron const& poly() const {return m_poly;}
    
    /// Grows a convex shape to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      m_poly.add_generator(point_generator(point));
    }
    
    /// Grows a convex shape to include the given convex shape.
    template <class RealT1>
    void grow( BConvex<RealT1> const& bconv ) {
      m_poly.poly_hull_assign(bconv.poly());
    }
    
    /// Crops (intersects) this convex shape to the given convex shape.
    template <class RealT1>
    void crop( BConvex<RealT1> const& bconv ) {
      m_poly.intersection_assign(bconv.poly());
    }
    
    /// Expands this convex shape by the given offset in every direction.
    void expand( RealT offset ) { //FIXME
      VW_ASSERT(0, LogicErr() << "expand() is not implemented for BConvex!");
    }

    /// Contracts this convex shape by the given offset in every direction.
    void contract( RealT offset ) { //FIXME
      VW_ASSERT(0, LogicErr() << "contract() is not implemented for BConvex!");
    }
    
    /// Returns true if the given point is contained in the convex shape.
    template <class VectorT>
    bool contains( VectorBase<VectorT> const& point ) const {
      return (m_poly.relation_with(point_generator(point)) == Parma_Polyhedra_Library::Poly_Gen_Relation::subsumes());
    }
    
    /// Returns true if the given convex shape is entirely contained
    /// in this convex shape.
    bool contains( BConvex const& bconv ) const {
      return m_poly.contains(bconv.poly());
    }
    
    /// Returns true if the given convex shape intersects this
    /// convex shape.
    bool intersects( BConvex const& bconv ) const {
      return !m_poly.is_disjoint_from(bconv.poly());
    }
    
    /// Returns the size (i.e. twice the largest distance between
    /// the center and a vertex) of the convex shape.
    RealT size() const {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      double result = 0;
      double d;
      Vector<double> v(dim);
      Vector<double> c = center_();
      const Generator_System &gs = m_poly.minimized_generators();
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
        const Generator &g = *i;
        VW_ASSERT(g.is_point(), LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
        for (unsigned j = dim - 1; 1; j--) {
          v(j) = (RealT)(g.coefficient(Variable(j)).get_d());
          if (j == 0)
            break;
        }
        d = norm_2_sqr(v - c);
        result = std::max(result, d);
      }
      result = 2*std::sqrt(result);
      return (RealT)result;
    }

    /// Returns the center point of the convex shape.
    Vector<RealT> center() const {
      Vector<double> c = center_();
      Vector<RealT> result(c.size());
      for (unsigned i = 0; i < c.size(); i++)
        result(i) = (RealT)c(i);
      return result;
    }

    /// Returns true if the convex shape is empty (i.e. degenerate).
    bool empty() const {
      return m_poly.is_empty();
    }

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator*=( ScalarT s ) {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      for (unsigned i = dim - 1; 1; i--) {
        Linear_Expression e;
        e += Variable(i) * (s * 10000);
        m_poly.affine_image(Variable(i), e, 10000);
        if (i == 0)
          break;
      }
      return *this;
    }

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator/=( ScalarT s ) {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      for (unsigned i = dim - 1; 1; i--) {
        Linear_Expression e;
        e += Variable(i) * 10000;
        m_poly.affine_image(Variable(i), e, s * 10000);
        if (i == 0)
          break;
      }
      return *this;
    }

    /// Offsets the convex shape by the given vector.
    template <class VectorT>
    BConvex& operator+=( VectorBase<VectorT> const& v ) {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      for (unsigned i = dim - 1; 1; i--) {
        Linear_Expression e;
        e += (Variable(i) * 10000) + (v.impl()(i) * 10000);
        m_poly.affine_image(Variable(i), e, 10000);
        if (i == 0)
          break;
      }
      return *this;
    }

    /// Offsets the convex shape by the negation of the given vector.
    template <class VectorT>
    BConvex& operator-=( VectorBase<VectorT> const& v ) {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      for (unsigned i = dim - 1; 1; i--) {
        Linear_Expression e;
        e += (Variable(i) * 10000) + ((-v.impl()(i)) * 10000);
        m_poly.affine_image(Variable(i), e, 10000);
        if (i == 0)
          break;
      }
      return *this;
    }

  protected:
    /// Returns the center point of the convex shape.
    Vector<double> center_() const {
      using namespace Parma_Polyhedra_Library;
      unsigned dim = m_poly.space_dimension();
      Vector<double> c(dim);
      unsigned n = 0;
      const Generator_System &gs = m_poly.minimized_generators();
      fill(c, 0);
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
        const Generator &g = *i;
        VW_ASSERT(g.is_point(), LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
        for (unsigned j = dim - 1; 1; j--) {
          c(j) += g.coefficient(Variable(j)).get_d();
          if (j == 0)
            break;
        }
        n++;
      }
      c /= n;
      return c;
    }
  
    /// Creates a PPL point.
    template <class VectorT>
    Parma_Polyhedra_Library::Generator point_generator( VectorBase<VectorT> const& point ) const { //FIXME: non-constant precision
      using namespace Parma_Polyhedra_Library;
      //using namespace ::Parma_Polyhedra_Library::IO_Operators; // for debugging
      Linear_Expression e;
      for (unsigned i = point.impl().size() - 1; 1; i--) {
        e += (point.impl()(i) * 10000) * Variable(i);
        if (i == 0)
          break;
      }
      Generator g = ::Parma_Polyhedra_Library::point(e, 10000);
      //std::cout << e << std::endl;
      return g;
    }
    
    Parma_Polyhedra_Library::C_Polyhedron m_poly;
  };
  
  /// Scales a convex shape relative to the origin.
  template <class RealT, class ScalarT>
  inline BConvex<RealT> operator*( BConvex<RealT> const& bconv, ScalarT s ) {
    BConvex<RealT> result = bconv;
    result *= s;
    return result;
  }

  /// Scales a convex shape relative to the origin.
  template <class RealT, class ScalarT>
  inline BConvex<RealT> operator/( BConvex<RealT> const& bconv, ScalarT s ) {
    BConvex<RealT> result = bconv;
    result /= s;
    return result;
  }

  /// Scales a convex shape relative to the origin.
  template <class RealT, class ScalarT>
  inline BConvex<RealT> operator*( ScalarT s, BConvex<RealT> const& bconv ) {
    return bconv * s;
  }
  
  /// Offsets a convex shape by the given vector.
  template <class RealT, class VectorT>
  inline BConvex<RealT> operator+( BConvex<RealT> const& bconv, VectorBase<VectorT> const& v ) {
    BConvex<RealT> result = bconv;
    result += v.impl();
    return result;
  }

  /// Offsets a convex shape by the given vector.
  template <class RealT, class VectorT>
  inline BConvex<RealT> operator+( VectorBase<VectorT> const& v, BConvex<RealT> const& bconv ) {
    return bconv + v;
  }

  /// Offsets a convex shape by the negation of the given vector.
  template <class RealT, class VectorT>
  inline BConvex<RealT> operator-( BConvex<RealT> const& bconv, VectorBase<VectorT> const& v ) {
    BConvex<RealT> result = bconv;
    result -= v.impl();
    return result;
  }

  /// Equality of two convex shapes.
  template <class RealT1, class RealT2>
  inline bool operator==( BConvex<RealT1> const& bconv1, BConvex<RealT2> const& bconv2 ) {
    return bconv1.poly() == bconv2.poly();
  }

  /// Inequality of two convex shapes.
  template <class RealT1, class RealT2>
  inline bool operator!=( BConvex<RealT1> const& bconv1, BConvex<RealT2> const& bconv2 ) {
    return bconv1.poly() != bconv2.poly();
  }
  
  /// Writes a convex shape to an ostream.
  template <class RealT>
  std::ostream& operator<<( std::ostream& os, BConvex<RealT> const& bconv ) {
    using ::Parma_Polyhedra_Library::IO_Operators::operator<<;
    return os << bconv.poly();
  }

} // namespace math

  // Convenience typedefs
  using math::BConvex;
  typedef BConvex<float64> BConvexN;
  typedef BConvex<float32> BConvexNf;
  typedef BConvex<int32> BConvexNi;
} // namespace vw

#endif // __VW_MATH__BCONVEX_H__
