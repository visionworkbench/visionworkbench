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

#include <iostream>
#include <math.h>
#include <vector>

#include <vw/Math/Vector.h>
#include <vw/Math/BConvex.h>
using namespace vw::math::bconvex_rational;

#include <ppl.hh> // Parma Polyhedra Library
using namespace Parma_Polyhedra_Library;

namespace {
  /// Creates a PPL point.
  Generator point_generator( vw::math::Vector<Rational> const& point ) {
    Linear_Expression e;
    if (point(0).is_signed) {
      for (unsigned i = point.size() - 1; 1; i--) {
        e += point(i).signed_num * Variable(i);
        if (i == 0)
          break;
      }
    }
    else {
      for (unsigned i = point.size() - 1; 1; i--) {
        e += point(i).unsigned_num * Variable(i);
        if (i == 0)
          break;
      }
    }
    Generator g = ::Parma_Polyhedra_Library::point(e, point(0).den);
    return g;
  }
  
  /// Multiplies a PPL Variable by a scalar.
  template <class ScalarT>
  Linear_Expression mult_variable(Variable v, ScalarT s) {
    return v * s;
  }
  
  /// Adds a scalar to a PPL Linear_Expression.
  template <class ScalarT>
  Linear_Expression add_scalar(Linear_Expression e, ScalarT s) {
    return e + s;
  }
  
  /// Subtracts a scalar from a PPL Linear_Expression.
  template <class ScalarT>
  Linear_Expression sub_scalar(Linear_Expression e, ScalarT s) {
    return e - s;
  }

} // namespace

namespace vw {
namespace math {

  /*static*/ void *BConvex::new_poly( unsigned dim ) {
    C_Polyhedron *p = new C_Polyhedron( dim, EMPTY );
    return (void*)p;
  }

  /*static*/ void *BConvex::new_poly( const void *poly ) {
    C_Polyhedron *p = new C_Polyhedron( *((const C_Polyhedron*)poly) );
    return (void*)p;
  }
  
  /*static*/ void BConvex::delete_poly( void *poly ) {
    if (poly)
      delete (C_Polyhedron*)poly;
  }
  
  /*static*/ void BConvex::copy_poly( const void *from_poly, void *to_poly ) {
    *((C_Polyhedron*)to_poly) = *((const C_Polyhedron*)from_poly);
  }

  void BConvex::grow_( Vector<Rational> const& point ) {
    ((C_Polyhedron*)m_poly)->add_generator(point_generator(point));
  }

  void BConvex::grow( BConvex const& bconv ) {
    ((C_Polyhedron*)m_poly)->poly_hull_assign(*((const C_Polyhedron*)(bconv.poly())));
  }

  void BConvex::crop( BConvex const& bconv ) {
    ((C_Polyhedron*)m_poly)->intersection_assign(*((const C_Polyhedron*)(bconv.poly())));
  }

  bool BConvex::contains_( Vector<Rational> const& point ) const {
    return (((const C_Polyhedron*)m_poly)->relation_with(point_generator(point)) == Poly_Gen_Relation::subsumes());
  }

  bool BConvex::contains( BConvex const& bconv ) const {
    return ((const C_Polyhedron*)m_poly)->contains(*((const C_Polyhedron*)(bconv.poly())));
  }

  bool BConvex::intersects( BConvex const& bconv ) const {
    return !((const C_Polyhedron*)m_poly)->is_disjoint_from(*((const C_Polyhedron*)(bconv.poly())));
  }

  double BConvex::size() const {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    double result = 0;
    double d;
    Vector<double> v(dim);
    Vector<double> c = center();
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      const Generator &g = *i;
      VW_ASSERT(g.is_point(), LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
      for (unsigned j = dim - 1; 1; j--) {
        v(j) = g.coefficient(Variable(j)).get_d();
        if (j == 0)
          break;
      }
      d = norm_2_sqr(v - c);
      result = std::max(result, d);
    }
    result = 2*std::sqrt(result);
    return result;
  }

  Vector<double> BConvex::center() const {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<double> c(dim);
    unsigned n = 0;
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
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

  bool BConvex::equal( BConvex const& bconv ) const {
    return *((const C_Polyhedron*)m_poly) == *((const C_Polyhedron*)(bconv.poly()));
  }

  bool BConvex::empty() const {
    return ((const C_Polyhedron*)m_poly)->is_empty();
  }

  void BConvex::print( std::ostream& os ) const {
    using ::Parma_Polyhedra_Library::IO_Operators::operator<<;
    os << *((const C_Polyhedron*)m_poly);
  }

  void BConvex::operator_mult_eq_( Rational const& s ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      if (s.is_signed)
        e += mult_variable(Variable(i), s.signed_num);
      else
        e += mult_variable(Variable(i), s.unsigned_num);
      ((C_Polyhedron*)m_poly)->affine_image(Variable(i), e, s.den);
      if (i == 0)
        break;
    }
  }

  void BConvex::operator_div_eq_( Rational const& s ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      e += mult_variable(Variable(i), s.den);
      if (s.is_signed)
        ((C_Polyhedron*)m_poly)->affine_image(Variable(i), e, s.signed_num);
      else
        ((C_Polyhedron*)m_poly)->affine_image(Variable(i), e, s.unsigned_num);
      if (i == 0)
        break;
    }
  }

  void BConvex::operator_plus_eq_( Vector<Rational> const& v ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      if (v(0).is_signed)
        e += add_scalar(mult_variable(Variable(i), v(i).den), v(i).signed_num);
      else
        e += add_scalar(mult_variable(Variable(i), v(i).den), v(i).unsigned_num);
      ((C_Polyhedron*)m_poly)->affine_image(Variable(i), e, v(0).den);
      if (i == 0)
        break;
    }
  }

  void BConvex::operator_minus_eq_( Vector<Rational> const& v ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      if (v(0).is_signed)
        e += sub_scalar(mult_variable(Variable(i), v(i).den), v(i).signed_num);
      else
        e += sub_scalar(mult_variable(Variable(i), v(i).den), v(i).unsigned_num);
      ((C_Polyhedron*)m_poly)->affine_image(Variable(i), e, v(0).den);
      if (i == 0)
        break;
    }
  }

}} // namespace vw::math
