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
#include <vw/Math/Matrix.h>
#include <vw/Math/BConvex.h>
using namespace vw::math::bconvex_rational;

#include <ppl.hh> // Parma Polyhedra Library
using namespace Parma_Polyhedra_Library;

#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
extern "C" {
#include <qhull/qhull.h>
#include <qhull/qset.h>
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

  /// Returns the PPL Generator (point) as a Vector.
  void convert_point( Generator const& g, vw::Vector<double> &p ) {
    VW_ASSERT(g.is_point(), vw::LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
    unsigned dim = g.space_dimension();
    unsigned i;
    p.set_size(dim);
    for (i = dim - 1; 1; i--) {
      p(i) = g.coefficient(Variable(i)).get_d() / g.divisor().get_d();
      if (i == 0)
        break;
    }
  }

  /// Finds the number of constraints in a constraint system.
  inline unsigned num_constraints( Constraint_System const& cs ) {
    unsigned n = 0;
    Constraint_System::const_iterator i;
    for (i = cs.begin(); i != cs.end(); i++)
      n++;
    return n;
  }

  /// Finds the number of generators in a generator system.
  inline unsigned num_generators( Generator_System const& gs ) {
    unsigned n = 0;
    Generator_System::const_iterator i;
    for (i = gs.begin(); i != gs.end(); i++)
      n++;
    return n;
  }

#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1

  /// Runs qhull.
  void qhull_run( unsigned dim, unsigned num_points, double *p ) {
    int retval;
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
  }

  /// Frees qhull memory.
  void qhull_free() {
    int curlong, totlong;
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
    VW_ASSERT(curlong == 0 && totlong == 0, vw::LogicErr() << "qhull did not free all of its memory");
  }

  /// Undo projection to plane.
  void write_vrml_unproject( unsigned dim, vertexT *vertex, vw::Matrix<double> const& basis, vw::Vector<double> const& plane_origin, vw::Vector<double> &p ) {
    unsigned i;
    if (dim == 2) {
      p.set_size(2);
      for (i = 0; i < dim; i++)
        p[i] = vertex->point[i];
    }
    else {
      p = plane_origin;
      for (i = 0; i < dim - 1; i++)
        p += elem_prod(vertex->point[i], select_col(basis, i + 1));
    }
  }

  /// Write a point in VRML.
  void write_vrml_point( std::ostream& os, vw::Vector<double> const& p ) {
    unsigned i;
    os << "   Separator {\n";
    os << "      Coordinate3 {\n";
    os << "         point [\n";
    os << "          ";
    for (i = 0; i < p.size(); i++)
      os << " " << p[i];
    for (; i < 3; i++)
      os << " 0";
    os << ",\n";
    os << "         ]\n";
    os << "      }\n";
    os << "      PointSet {\n";
    os << "         startIndex 0\n";
    os << "         numPoints 1\n";
    os << "      }\n";
    os << "   }\n";
  }

  /// Write a line in VRML.
  void write_vrml_line( std::ostream& os, vw::Vector<double> const& p1, vw::Vector<double> const& p2 ) {
    unsigned i;
    os << "   Separator {\n";
    os << "      Coordinate3 {\n";
    os << "         point [\n";
    os << "          ";
    for (i = 0; i < p1.size(); i++)
      os << " " << p1[i];
    for (; i < 3; i++)
      os << " 0";
    os << ",\n";
    os << "          ";
    for (i = 0; i < p2.size(); i++)
      os << " " << p2[i];
    for (; i < 3; i++)
      os << " 0";
    os << ",\n";
    os << "         ]\n";
    os << "      }\n";
    os << "      IndexedLineSet {\n";
    os << "         coordIndex [ 0, 1 ]\n";
    os << "      }\n";
    os << "   }\n";
  }

#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1

  /// Find an orthonormal basis using a given Vector.
  void gram_schmidt( vw::Vector<double> const& v, vw::Matrix<double> &b ) {
    using namespace vw;
    using namespace ::vw::math;
    unsigned dim = v.size();
    unsigned largest;
    unsigned row, col, col2;
    double d;
    b.set_size(dim, dim);
    select_col(b, 0) = v;
    normalize(select_col(b, 0));
    largest = index_norm_inf(select_col(b, 0));
    for (row = 0, col = 1; col < dim; row++, col++) {
      if (row == largest)
        row++;
      select_col(b, col) = select_col(b, 0);
      b(row, col) += 1.0;
      for (col2 = 0; col2 < col; col2++) {
        d = dot_prod(select_col(b, col), select_col(b, col2));
        select_col(b, col) -= elem_prod(d, select_col(b, col2));
      }
      normalize(select_col(b, col));
      //for (col2 = 0; col2 < col; col2++) {
      //  VW_ASSERT(dot_prod(select_col(b, col), select_col(b, col2)) == 0, LogicErr() << "gram_schmidt failed!");
      //}
    }
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
    facetT *facet;
    Vector<double> facetv(dim + 1);
    unsigned i;

    qhull_run(dim, num_points, p);
    
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
    
    qhull_free();
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
      coefficients(*i, coef);
      norm = norm_2(subvector(coef, 0, dim));
      coef /= norm;
      VW_ASSERT(coef[dim] >= 0, LogicErr() << "Center is outside of polyhedron!");
      coef[dim] += offset;
      Constraint c = constraint_generator(coef);
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
      convert_point(*i, v);
      d = norm_2_sqr(v - c);
      result = std::max(result, d);
    }
    result = 2*std::sqrt(result);
    return result;
  }

  Vector<double> BConvex::center() const {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<double> c(dim), v(dim);
    unsigned n = 0;
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    fill(c, 0);
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, v);
      c += v;
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

  void BConvex::write_vrml( std::ostream& os ) const {
#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    unsigned num_points = num_generators(gs);
    Vector<double> vertexv( dim ), vertexv2( dim );

    VW_ASSERT(dim > 0 && dim <= 3, ArgumentErr() << "Cannot write vrml of dimension <= 0 or > 3!");

    os << "#VRML V1.0 ascii\n\n";
    os << "# Created by the Intelligent Robotics Group,\n";
    os << "# NASA Ames Research Center\n";
    os << "# File generated by the NASA Ames Vision Workbench.\n\n";
    os << "Separator {\n";
    os << "   Material {\n";
    os << "      ambientColor 1.00 1.00 1.00\n";
    os << "      diffuseColor 1.00 1.00 1.00\n";
    os << "      specularColor 0.00 0.00 0.00\n";
    os << "      emissiveColor 0.00 0.00 0.00\n";
    os << "      shininess 0.00\n";
    os << "      transparency 0.00\n";
    os << "   }\n";

    if (dim == 1) {
      VW_ASSERT(num_points <= 2, LogicErr() << "Found line with more than 2 vertices!");
      if (num_points == 1) {
        convert_point(*(gs.begin()), vertexv);
        write_vrml_point(os, vertexv);
      }
      else {
        convert_point(*(gs.begin()), vertexv);
        convert_point(*(++(gs.begin())), vertexv2);
        write_vrml_line(os, vertexv, vertexv2);
      }
    }
    else if (dim == 2 || dim == 3){
      double *p = new double[dim * num_points];
      double **ps = 0;
      facetT *facet;
      Vector<double> facetv( dim + 1 );
      unsigned num_facets;
      vertexT *vertex;
      vertexT **vertexp;
      unsigned *num_vertices = 0;
      unsigned num_vertices_;
      Matrix<double> *basis = 0;
      Vector<double> *plane_origin = 0;
      Vector<double> proj;
      unsigned j, k, l;
  
      k = 0;
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
        convert_point(*i, vertexv);
        for (unsigned j = 0; j < dim; j++, k++)
          p[k] = vertexv[j];
      }
      VW_ASSERT(k == dim * num_points, LogicErr() << "k != dim * num_points!");
  
      qhull_run(dim, num_points, p);
      if (dim == 2) {
        num_facets = 1;
        ps = new double*[num_facets];
        num_vertices = new unsigned[num_facets];
        basis = new Matrix<double>[num_facets];
        plane_origin = new Vector<double>[num_facets];
        j = 0;
        num_vertices[j] = 0;
        FORALLvertices {
          num_vertices[j]++;
        }
        ps[j] = new double[dim * num_vertices[j]];
        k = 0;
        FORALLvertices {
          for (l = 0; l < dim; l++)
            vertexv[l] = vertex->point[l];
          for (l = 0; l < dim; l++, k++)
            ps[j][k] = proj[l];
        }
      }
      else {
        num_facets = 0;
        FORALLfacets {
          num_facets++;
        }
        ps = new double*[num_facets];
        num_vertices = new unsigned[num_facets];
        basis = new Matrix<double>[num_facets];
        plane_origin = new Vector<double>[num_facets];
        j = 0;
        FORALLfacets {
          num_vertices[j] = 0;
          FOREACHvertex_(facet->vertices) {
            num_vertices[j]++;
          }
          ps[j] = new double[(dim - 1) * num_vertices[j]];
          for (l = 0; l < dim; l++)
            facetv[l] = -facet->normal[l];
          facetv[l] = -facet->offset;
          gram_schmidt(subvector(facetv, 0, dim), basis[j]);
          facetv /= norm_2(subvector(facetv, 0, dim));
          plane_origin[j] = elem_prod(-facetv[dim], subvector(facetv, 0, dim));
          k = 0;
          FOREACHvertex_(facet->vertices) {
            for (l = 0; l < dim; l++)
              vertexv[l] = vertex->point[l];
            proj = transpose(basis[j]) * vertexv;
            //VW_ASSERT(proj[0] == -facetv[dim], LogicErr() << "Projection to plane failed!");
            for (l = 1; l < dim; l++, k++)
              ps[j][k] = proj[l];
          }
          j++;
        }
      }
      qhull_free();
  
      for (j = 0; j < num_facets; j++) {
        qhull_run(2, num_vertices[j], ps[j]);
        FORALLfacets {
          num_vertices_ = 0;
          FOREACHvertex_(facet->vertices) {
            num_vertices_++;
          }
          VW_ASSERT(num_vertices_ <= 2, LogicErr() << "Found facet with more than 2 vertices!");
          if (num_vertices_ == 1) {
            vertex = SETfirstt_(facet->vertices, vertexT);
            write_vrml_unproject(dim, vertex, basis[j], plane_origin[j], vertexv);
            write_vrml_point(os, vertexv);
          }
          else if (num_vertices_ == 2) {
            vertex = SETfirstt_(facet->vertices, vertexT);
            write_vrml_unproject(dim, vertex, basis[j], plane_origin[j], vertexv);
            vertex = SETsecondt_(facet->vertices, vertexT);
            write_vrml_unproject(dim, vertex, basis[j], plane_origin[j], vertexv2);
            write_vrml_line(os, vertexv, vertexv2);
          }
        }
        qhull_free();
      }
  
      delete[] p;
      p = 0;
      for (j = 0; j < num_facets; j++)
        delete[] ps[j];
      delete[] ps;
      ps = 0;
      delete[] num_vertices;
      num_vertices = 0;
      delete[] basis;
      basis = 0;
      delete[] plane_origin;
      plane_origin = 0;
    }

    os << "}\n";

#else // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    VW_ASSERT(0, ArgumentErr() << "Must have qhull to print vrml!");
#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
  }
  
  unsigned BConvex::num_facets() const {
    return num_constraints(((const C_Polyhedron*)m_poly)->minimized_constraints());
  }
  
  BBoxN BConvex::bounding_box() const {
    BBoxN bbox;
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<double> v(dim);
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, v);
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
