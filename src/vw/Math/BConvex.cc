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

#include <boost/thread/mutex.hpp>

#include <vw/config.h> // VW_HAVE_PKG_QHULL
#include <vw/Math/Vector.h>
#include <vw/Math/BConvex.h>
using namespace vw::math::bconvex_rational;

#include <ppl.hh> // Parma Polyhedra Library
using namespace Parma_Polyhedra_Library;

#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
extern "C" {
#include <qhull/qhull.h>
}
#endif

namespace {
  /// Mutex for qhull, which is not thread-safe.
  boost::mutex bconvex_qhull_mutex;

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
  inline Linear_Expression mult_variable(Variable v, ScalarT s) {
    return v * s;
  }
  
  /// Adds a scalar to a PPL Linear_Expression.
  template <class ScalarT>
  inline Linear_Expression add_scalar(Linear_Expression e, ScalarT s) {
    return e + s;
  }
  
  /// Subtracts a scalar from a PPL Linear_Expression.
  template <class ScalarT>
  inline Linear_Expression sub_scalar(Linear_Expression e, ScalarT s) {
    return e - s;
  }

  /// Creates a PPL Constraint from a Vector.
  Constraint constraint_generator( vw::math::Vector<double> const& coef ) {
    unsigned dim = coef.size() - 1;
    vw::math::Vector<Rational> coefr( coef.size() );
    Linear_Expression e;
    unsigned i;
    vw::math::BConvex::convert_vector(coef, coefr);
    for (i = dim - 1; 1; i--) {
      if (coefr[0].is_signed)
        e += mult_variable(Variable(i), coefr[i].signed_num);
      else
        e += mult_variable(Variable(i), coefr[i].unsigned_num);
      if (i == 0)
        break;
    }
    if (coefr[0].is_signed)
      e += coefr[dim].signed_num;
    else
      e += coefr[dim].unsigned_num;
    Constraint c = (e >= 0);
    return c;
  }

  /// Finds the coefficients of a PPL Constraint.
  void coefficients( Constraint const& c, vw::math::Vector<double> &coef ) {
    unsigned dim = c.space_dimension();
    unsigned i;
    coef.set_size(dim + 1);
    for (i = dim - 1; 1; i--) {
      coef[i] = c.coefficient(Variable(i)).get_d();
      if (i == 0)
        break;
    }
    coef[dim] = c.inhomogeneous_term().get_d();
  }

} // namespace

namespace vw {
namespace math {

  namespace bconvex_rational {

    std::ostream& operator<<( std::ostream& os, Rational const& r ) {
      if (r.is_signed)
        return os << r.signed_num << "/" << r.den;
      return os << r.unsigned_num << "/" << r.den;
    }

  } // namespace bconvex_rational

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
  
  /*static*/ bool BConvex::have_qhull() {
#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    return true;
#else
    return false;
#endif
  }

  void BConvex::init_with_qhull( unsigned dim, unsigned num_points, double *p ) {
#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    boost::mutex::scoped_lock lock(bconvex_qhull_mutex);
    int retval;
    int curlong, totlong;
    facetT *facet;
    Vector<double> facetv(dim + 1);
    FILE *fake_stdout;
    FILE *fake_stderr;
    unsigned i;
  
    fake_stdout = tmpfile();
    if (!fake_stdout)
      fake_stdout = stdout;
    fake_stderr = tmpfile();
    if (!fake_stderr)
      fake_stderr = stderr;
    retval = qh_new_qhull((int)dim, (int)num_points,
                          p, False, "qhull s Tcv", fake_stdout, fake_stderr);
    if (retval != 0 && fake_stderr != stderr) {
      std::string qhull_error;
      int c;
      rewind(fake_stderr);
      while ((c = fgetc(fake_stderr)) != EOF)
        qhull_error.append(1, (char)c);
      VW_ASSERT(retval == 0, vw::LogicErr() << "qhull returned with error: " << qhull_error);
    }
    VW_ASSERT(retval == 0, vw::LogicErr() << "qhull returned with error.");
    if (fake_stdout != stdout)
      fclose(fake_stdout);
    if (fake_stderr != stderr)
      fclose(fake_stderr);
    
    m_poly = (void*)(new C_Polyhedron( dim, UNIVERSE ));
    FORALLfacets {
      for (i = 0; i < dim; i++)
        facetv[i] = -facet->normal[i];
      facetv[i] = -facet->offset;
      //NOTE: Printing facet->toporient here for the [0,1]^3 unit cube test case (see TestBConvex.h) demonstrates that facet->toporient is meaningless in the qhull output.
      Constraint c = constraint_generator(facetv);
      ((C_Polyhedron*)m_poly)->add_constraint(c);
    }
    
    VW_ASSERT(!this->empty(), LogicErr() << "Convex hull is empty!");
    VW_ASSERT(!((const C_Polyhedron*)m_poly)->is_universe(), LogicErr() << "Convex hull is universe!");
    
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
    VW_ASSERT(curlong == 0 && totlong == 0, vw::LogicErr() << "qhull did not free all of its memory");
#else // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    return;
#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
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

  void BConvex::expand( double offset ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<double> center_ = center();
    C_Polyhedron *poly = new C_Polyhedron( dim, UNIVERSE );
    Vector<double> coef( dim + 1 );
    double norm;
    *this -= center_;
    const Constraint_System &cs = ((const C_Polyhedron*)m_poly)->minimized_constraints();
    Constraint_System::const_iterator i;
    for (i = cs.begin(); i != cs.end(); i++) {
      {
        using ::Parma_Polyhedra_Library::IO_Operators::operator<<;
        ::std::cout << "Original: " << *i << std::endl;
      }
      coefficients(*i, coef);
      norm = norm_2(subvector(coef, 0, dim));
      coef /= norm;
      VW_ASSERT(coef[dim] >= 0, LogicErr() << "Center is outside of polyhedron!");
      coef[dim] += offset;
      Constraint c = constraint_generator(coef);
      {
        using ::Parma_Polyhedra_Library::IO_Operators::operator<<;
        ::std::cout << "Modified: " << c << std::endl;
      }
      poly->add_constraint(c);
    }
    delete (C_Polyhedron*)m_poly;
    m_poly = (void*)poly;
    *this += center_;
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
        v(j) = g.coefficient(Variable(j)).get_d() / g.divisor().get_d();
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
        c(j) += g.coefficient(Variable(j)).get_d() / g.divisor().get_d();
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
  
  unsigned BConvex::num_facets() const {
    unsigned n = 0;
    const Constraint_System &cs = ((const C_Polyhedron*)m_poly)->minimized_constraints();
    Constraint_System::const_iterator i;
    for (i = cs.begin(); i != cs.end(); i++)
      n++;
    return n;
  }
  
  BBoxN BConvex::bounding_box() const {
    BBoxN bbox;
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<double> v(dim);
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      const Generator &g = *i;
      VW_ASSERT(g.is_point(), LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
      for (unsigned j = dim - 1; 1; j--) {
        v(j) = g.coefficient(Variable(j)).get_d() / g.divisor().get_d();
        if (j == 0)
          break;
      }
      bbox.grow(v);
    }
    return bbox;
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

  std::ostream& operator<<( std::ostream& os, BConvex const& bconv ) {
    bconv.print(os);
    return os;
  }

}} // namespace vw::math
