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
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

#include <boost/thread/mutex.hpp>

#include <vw/config.h> // VW_HAVE_PKG_QHULL
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/PointListIO.h>
#include <vw/Math/SpatialTree.h>
#include <vw/Math/BConvex.h>
using namespace vw::math::bconvex_promote;

#include <gmpxx.h> // GNU Multiple Precision Arithmetic Library
#include <ppl.hh> // Parma Polyhedra Library
using namespace Parma_Polyhedra_Library;

#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
extern "C" {
#include <qhull/qhull.h>
#include <qhull/qset.h>
}
#endif

namespace vw {
namespace math {

  /// Square of vector 2-norm (with gmp elements).
  template <class VectorT>
  inline double norm_2_sqr_gmp( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) result += (*i).get_d() * (*i).get_d();
    return result;
  }

  /// Vector 2-norm (with gmp elements).
  template <class VectorT>
  inline double norm_2_gmp( VectorBase<VectorT> const& v ) {
    return sqrt( norm_2_sqr_gmp(v) );
  }

}} // namespace vw::math

namespace {
  /// Mutex for qhull, which is not thread-safe.
  boost::mutex bconvex_qhull_mutex;
  
  /// Point primitive for SpatialTree.
  class PointPrimitive : public vw::math::GeomPrimitive
  {
  public:
    PointPrimitive() {}
    
    ~PointPrimitive() {}
  
    // Implementation of GeomPrimitive interface.
    virtual double distance(const vw::Vector<double> &point_) const
    {
      using namespace vw::math;
      return norm_2(point - point_);
    }
    virtual const vw::BBox<double> &bounding_box() const {return bbox;}
  
    unsigned index;
    vw::Vector<double> point;
    vw::BBox<double> bbox;
  };
  
  /// Indexed polyhedron edge.
  struct IndexedEdge {
    unsigned index[2];
    bool used;
  };

  /// Find the least common denominator.
  mpz_class least_common_denominator( vw::math::Vector<mpq_class> const& v ) {
    unsigned dim = v.size();
    VW_ASSERT(dim > 0, vw::LogicErr() << "Cannot find least common denominator of 0 points!");
    mpz_class den = v(0).get_den();
    for (unsigned i = 1; i < dim; i++) {
      mpq_class q(v(i).get_den(), den);
      q.canonicalize();
      den = q.get_den() * v(i).get_den();
    }
    return den;
  }

  /// Convert all rationals in Vector so that they have a common denominator.
  //NOTE: The resulting rationals may not be canonical, and arithmetic expressions
  //  involving them may not work properly. This is desired behavior, as we want
  //  all rationals in the Vector to have the same denominator.
  vw::math::Vector<mpq_class> const& common_denominator( vw::math::Vector<mpq_class> const& point, vw::math::Vector<mpq_class> &p ) {
    unsigned dim = point.size();
    p.set_size(dim);
    if (dim > 0) {
      mpz_class lcd = least_common_denominator(point);
      for (unsigned i = 0; i < dim; i++) {
        p(i).get_num() = point(i).get_num() * (lcd / point(i).get_den());
        p(i).get_den() = lcd;
      }
    }
    return p;
  }

  /// Convert all rationals in Vector so that they have no denominator
  /// (so that linear expression with rationals as coefficients is equivalent).
  vw::math::Vector<mpz_class> const& no_denominator( vw::math::Vector<mpq_class> const& point, vw::math::Vector<mpz_class> &p ) {
    unsigned dim = point.size();
    p.set_size(dim);
    if (dim > 0) {
      mpz_class lcd = least_common_denominator(point);
      for (unsigned i = 0; i < dim; i++) {
        p(i) = point(i).get_num() * (lcd / point(i).get_den());
      }
    }
    return p;
  }

  /// Converts a Promoted Vector to an mpq_class Vector.
  vw::math::Vector<mpq_class> const& convert_vector( vw::math::Vector<Promoted> const& point, vw::math::Vector<mpq_class> &p ) {
    unsigned dim = point.size();
    p.set_size(dim);
    if (!point(0).is_integral) {
      for (unsigned i = 0; i < dim; i++)
        p(i) = mpq_class(point(i).val.f);
    }
    else if (point(0).is_signed) {
      for (unsigned i = 0; i < dim; i++)
        p(i) = mpq_class(point(i).val.s);
    }
    else {
      for (unsigned i = 0; i < dim; i++)
        p(i) = mpq_class(point(i).val.u);
    }
    return p;
  }

  /// Converts a double Vector to an mpq_class Vector.
  vw::math::Vector<mpq_class> const& convert_vector( vw::math::Vector<double> const& point, vw::math::Vector<mpq_class> &p ) {
    unsigned dim = point.size();
    p.set_size(dim);
    for (unsigned i = 0; i < dim; i++)
      p(i) = mpq_class(point(i));
    return p;
  }

  /// Converts an mpq_class Vector to a double Vector.
  vw::math::Vector<double> const& unconvert_vector( vw::math::Vector<mpq_class> const& point, vw::math::Vector<double> &p ) {
    unsigned dim = point.size();
    p.set_size(dim);
    for (unsigned i = 0; i < dim; i++)
      p(i) = point(i).get_d();
    return p;
  }

  /// Converts a Promoted to an mpq_class.
  mpq_class const& convert_scalar( Promoted const& scalar, mpq_class &s ) {
    if (!scalar.is_integral) {
      s = mpq_class(scalar.val.f);
    }
    else if (scalar.is_signed) {
      s = mpq_class(scalar.val.s);
    }
    else {
      s = mpq_class(scalar.val.u);
    }
    return s;
  }

  /// Converts a double to an mpq_class.
  mpq_class const& convert_scalar( double const& scalar, mpq_class &s ) {
    s = mpq_class(scalar);
    return s;
  }

  /// Converts an mpq_class to a double.
  double const& unconvert_scalar( mpq_class const& scalar, double &s ) {
    s = scalar.get_d();
    return s;
  }

  /// Creates a PPL point.
  Generator point_generator( vw::math::Vector<mpq_class> const& point ) {
    Linear_Expression e;
    unsigned dim = point.size();
    vw::math::Vector<mpq_class> pointc( dim );
    common_denominator(point, pointc);
    for (unsigned i = point.size() - 1; 1; i--) {
      e += pointc(i).get_num() * Variable(i);
      if (i == 0)
        break;
    }
    Generator g = ::Parma_Polyhedra_Library::point(e, pointc(0).get_den());
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

  /// Creates a PPL Constraint from an mpz_class Vector.
  Constraint constraint_generator( vw::math::Vector<mpz_class> const& coef ) {
    unsigned dim = coef.size() - 1;
    Linear_Expression e;
    unsigned i;
    for (i = dim - 1; 1; i--) {
      e += mult_variable(Variable(i), coef[i]);
      if (i == 0)
        break;
    }
    e += coef[dim];
    Constraint c = (e >= 0);
    return c;
  }

  /// Creates a PPL Constraint from an mpq_class Vector.
  Constraint constraint_generator( vw::math::Vector<mpq_class> const& coef ) {
    vw::math::Vector<mpz_class> coefz( coef.size() );
    no_denominator(coef, coefz);
    return constraint_generator(coefz);
  }

  /// Finds the coefficients of a PPL Constraint.
  void coefficients( Constraint const& c, vw::math::Vector<mpz_class> &coef ) {
    unsigned dim = c.space_dimension();
    unsigned i;
    coef.set_size(dim + 1);
    for (i = dim - 1; 1; i--) {
      coef[i] = c.coefficient(Variable(i));
      if (i == 0)
        break;
    }
    coef[dim] = c.inhomogeneous_term();
  }

  /// Finds the coefficients of a PPL Constraint.
  void coefficients( Constraint const& c, vw::math::Vector<mpq_class> &coef ) {
    unsigned dim = c.space_dimension();
    unsigned i;
    coef.set_size(dim + 1);
    for (i = dim - 1; 1; i--) {
      coef[i] = mpq_class(c.coefficient(Variable(i)), 1);
      if (i == 0)
        break;
    }
    coef[dim] = mpq_class(c.inhomogeneous_term(), 1);
  }

  /// Returns the PPL Generator (point) as a Vector.
  void convert_point( Generator const& g, vw::Vector<mpq_class> &p ) {
    VW_ASSERT(g.is_point(), vw::LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
    unsigned dim = g.space_dimension();
    unsigned i;
    p.set_size(dim);
    for (i = dim - 1; 1; i--) {
      p(i) = mpq_class(g.coefficient(Variable(i)), g.divisor());
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
    for (i = gs.begin(); i != gs.end(); i++) {
      VW_ASSERT((*i).is_point(), vw::LogicErr() << "Retrieved non-point generator from closed convex polyhedron!");
      n++;
    }
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
  
#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1

  /// Find the center of a polyhedron.
  void poly_center( const C_Polyhedron *poly, vw::Vector<mpq_class> &c ) {
    VW_ASSERT(poly, vw::LogicErr() << "Polyhedron pointer is null!");
    unsigned dim = poly->space_dimension();
    vw::Vector<mpq_class> v(dim);
    unsigned n = 0;
    const Generator_System &gs = poly->minimized_generators();
    c.set_size(dim);
    fill(c, 0);
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, v);
      c += v;
      n++;
    }
    c /= n;
  }

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
    select_col(b, 0) /= norm_2(select_col(b, 0));
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
      select_col(b, col) /= norm_2(select_col(b, col));
      //for (col2 = 0; col2 < col; col2++) {
      //  VW_ASSERT(dot_prod(select_col(b, col), select_col(b, col2)) == 0, LogicErr() << "gram_schmidt failed!");
      //}
    }
  }

  /// Undo projection to plane.
  void plane_unproject( unsigned dim, vw::Vector<double> const& v, vw::Matrix<double> const& basis, vw::Vector<double> const& plane_origin, vw::Vector<double> &p ) {
    unsigned i;
    if (dim == 2) {
      p.set_size(2);
      for (i = 0; i < dim; i++)
        p[i] = v[i];
    }
    else {
      p = plane_origin;
      for (i = 0; i < dim - 1; i++)
        p += elem_prod(v[i], select_col(basis, i + 1));
    }
  }

  /// Find the facets of the polyhedron.
  void poly_facet_list( const C_Polyhedron *poly, std::vector<vw::Vector<mpq_class> > &points, std::vector<std::vector<unsigned> > &facets ) {
    if (!poly || poly->is_empty())
      return;
      
    unsigned dim = poly->space_dimension();
    const Generator_System &gs = poly->minimized_generators();
    unsigned num_points = num_generators(gs);
    vw::Vector<double> vertexv( dim );
    vw::Vector<mpq_class> vertexq( dim );

    VW_ASSERT(dim > 0 && dim <= 3, vw::ArgumentErr() << "Cannot find facet list of dimension <= 0 or > 3!");

    if (dim == 1) {
      VW_ASSERT(num_points <= 2, vw::LogicErr() << "Found line with more than 2 vertices!");
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
        convert_point(*i, vertexq);
        points.push_back(vertexq);
      }
      facets.push_back(std::vector<unsigned>());
      if (num_points == 1)
        facets[0].push_back(0);
      else {
        facets[0].push_back(0);
        facets[0].push_back(1);
      }
    }
    else if (dim == 2 || dim == 3) {
    
#if defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1

      boost::mutex::scoped_lock lock(bconvex_qhull_mutex);
      double *p = new double[dim * num_points];
      double **ps = 0;
      facetT *facet;
      vw::Vector<double> facetv( dim + 1 );
      unsigned num_facets;
      vertexT *vertex;
      vertexT **vertexp;
      vw::Vector<double> vertexv2( 2 );
      unsigned *num_vertices = 0;
      unsigned num_vertices_;
      vw::Matrix<double> *basis = 0;
      vw::Vector<double> *plane_origin = 0;
      vw::Vector<double> proj;
      vw::math::GeomPrimitive **prims;
      PointPrimitive *prim;
      IndexedEdge edge;
      std::vector<IndexedEdge> edges;
      unsigned edges_used;
      unsigned index;
      unsigned j, k, l;
      bool found_it;
  
      prims = new vw::math::GeomPrimitive*[num_points];
      k = 0;
      l = 0;
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++, l++) {
        prim = new PointPrimitive;
        prim->index = l;
        convert_point(*i, vertexq);
        points.push_back(vertexq);
        unconvert_vector(vertexq, prim->point);
        prim->bbox.grow(prim->point);
        for (j = 0; j < dim; j++, k++)
          p[k] = prim->point[j];
        prims[l] = (vw::math::GeomPrimitive*)prim;
      }
      VW_ASSERT(k == dim * num_points, vw::LogicErr() << "k != dim * num_points!");
      vw::math::SpatialTree st( num_points, prims, 5 );
  
      qhull_run(dim, num_points, p);
      if (dim == 2) {
        num_facets = 1;
        ps = new double*[num_facets];
        num_vertices = new unsigned[num_facets];
        basis = new vw::Matrix<double>[num_facets];
        plane_origin = new vw::Vector<double>[num_facets];
        j = 0;
        num_vertices[j] = 0;
        FORALLvertices {
          num_vertices[j]++;
        }
        ps[j] = new double[dim * num_vertices[j]];
        k = 0;
        FORALLvertices {
          for (l = 0; l < dim; l++, k++)
            ps[j][k] = vertex->point[l];
        }
      }
      else {
        num_facets = 0;
        FORALLfacets {
          num_facets++;
        }
        ps = new double*[num_facets];
        num_vertices = new unsigned[num_facets];
        basis = new vw::Matrix<double>[num_facets];
        plane_origin = new vw::Vector<double>[num_facets];
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
            //VW_ASSERT(proj[0] == -facetv[dim], vw::LogicErr() << "Projection to plane failed!");
            for (l = 1; l < dim; l++, k++)
              ps[j][k] = proj[l];
          }
          j++;
        }
      }
      qhull_free();
  
      for (j = 0; j < num_facets; j++) {
        facets.push_back(std::vector<unsigned>());
        qhull_run(2, num_vertices[j], ps[j]);
        edges.clear();
        FORALLfacets {
          num_vertices_ = 0;
          FOREACHvertex_(facet->vertices) {
            num_vertices_++;
          }
          VW_ASSERT(num_vertices_ <= 2, vw::LogicErr() << "Found facet with more than 2 vertices!");
          if (num_vertices_ == 1) {
            vertex = SETfirstt_(facet->vertices, vertexT);
            for (l = 0; l < 2; l++)
              vertexv2[l] = vertex->point[l];
            plane_unproject(dim, vertexv2, basis[j], plane_origin[j], vertexv);
            prim = (PointPrimitive*)st.closest(vertexv);
            VW_ASSERT(prim, vw::LogicErr() << "No closest vertex found!");
            VW_ASSERT(facets[j].empty(), vw::LogicErr() << "Found more than one singleton vertex!");
            facets[j].push_back(prim->index);
          }
          else if (num_vertices_ == 2) {
            vertex = SETfirstt_(facet->vertices, vertexT);
            for (l = 0; l < 2; l++)
              vertexv2[l] = vertex->point[l];
            plane_unproject(dim, vertexv2, basis[j], plane_origin[j], vertexv);
            prim = (PointPrimitive*)st.closest(vertexv);
            VW_ASSERT(prim, vw::LogicErr() << "No closest vertex found!");
            edge.index[0] = prim->index;
            vertex = SETsecondt_(facet->vertices, vertexT);
            for (l = 0; l < 2; l++)
              vertexv2[l] = vertex->point[l];
            plane_unproject(dim, vertexv2, basis[j], plane_origin[j], vertexv);
            prim = (PointPrimitive*)st.closest(vertexv);
            VW_ASSERT(prim, vw::LogicErr() << "No closest vertex found!");
            edge.index[1] = prim->index;
            edge.used = false;
            edges.push_back(edge);
          }
        }
        if (edges.empty()) {
          VW_ASSERT(!facets[j].empty(), vw::LogicErr() << "Found neither a singleton vertex nor an edge!");
        }
        else {
          VW_ASSERT(facets[j].empty(), vw::LogicErr() << "Found both a singleton vertex and an edge!");
          facets[j].push_back(edges[0].index[0]);
          index = edges[0].index[1];
          edges[0].used = true;
          edges_used = 1;
          while (edges_used < edges.size()) {
            facets[j].push_back(index);
            found_it = false;
            for (k = 1; k < edges.size(); k++) {
              if (!edges[k].used) {
                if (edges[k].index[0] == index) {
                  index = edges[k].index[1];
                  found_it = true;
                }
                else if (edges[k].index[1] == index) {
                  index = edges[k].index[0];
                  found_it = true;
                }
                if (found_it) {
                  edges[k].used = true;
                  edges_used++;
                  break;
                }
              }
            }
            VW_ASSERT(found_it, vw::LogicErr() << "Polygon is not closed!");
          }
          VW_ASSERT(index == edges[0].index[0], vw::LogicErr() << "Polygon is not closed!");
        }
        qhull_free();
      }
  
      delete[] p;
      p = 0;
      for (j = 0; j < num_facets; j++)
        delete[] ps[j];
      delete[] ps;
      ps = 0;
      for (j = 0; j < num_points; j++)
        delete (PointPrimitive*)prims[j];
      delete[] prims;
      prims = 0;
      delete[] num_vertices;
      num_vertices = 0;
      delete[] basis;
      basis = 0;
      delete[] plane_origin;
      plane_origin = 0;
      
#else // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
      VW_ASSERT(0, ArgumentErr() << "Must have qhull to find facet list for dimension >= 2!");
#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1

    }
  }

  /// Offsets the polyhedron by the given vector.
  void offset_poly( vw::Vector<mpq_class> const& v, C_Polyhedron *poly ) {
    unsigned dim = poly->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      e += add_scalar(mult_variable(Variable(i), v(i).get_den()), v(i).get_num());
      poly->affine_image(Variable(i), e, v(i).get_den());
      if (i == 0)
        break;
    }
  }

  /// Scales the polyhedron relative to the origin.
  void scale_poly( mpq_class const& s, C_Polyhedron *poly ) {
    unsigned dim = poly->space_dimension();
    for (unsigned i = dim - 1; 1; i--) {
      Linear_Expression e;
      e += mult_variable(Variable(i), s.get_num());
      poly->affine_image(Variable(i), e, s.get_den());
      if (i == 0)
        break;
    }
  }

} // namespace

namespace vw {
namespace math {

  namespace bconvex_promote {

    std::ostream& operator<<( std::ostream& os, Promoted const& r ) {
      if (!r.is_integral)
        return os << r.val.f;
      if (r.is_signed)
        return os << r.val.s;
      return os << r.val.u;
    }

  } // namespace bconvex_promote

  /*static*/ void *BConvex::new_poly( unsigned dim ) {
    C_Polyhedron *p = new C_Polyhedron( dim, EMPTY );
    return (void*)p;
  }

  /*static*/ void *BConvex::new_poly( const void *poly ) {
    C_Polyhedron *p = new C_Polyhedron( *((const C_Polyhedron*)poly) );
    return (void*)p;
  }

  /*static*/ void BConvex::new_poly_if_needed( unsigned dim, void *&poly ) {
    if (!poly)
      poly = new_poly(dim);
  }
  
  /*static*/ void BConvex::delete_poly( void *poly ) {
    if (poly)
      delete (C_Polyhedron*)poly;
  }
  
  /*static*/ void BConvex::copy_poly( const void *from_poly, void *&to_poly ) {
    if (to_poly)
      *((C_Polyhedron*)to_poly) = *((const C_Polyhedron*)from_poly);
    else
      to_poly = new_poly(from_poly);
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
    Vector<mpq_class> facetq(dim + 1);
    unsigned i;
    unsigned num_facets = 0, num_vertices = 0;
    vertexT *vertex;
    Vector<double> vertexv(dim);
    Vector<mpq_class> vertexq(dim);

    qhull_run(dim, num_points, p);
    
#if 0
    m_poly = (void*)(new C_Polyhedron( dim, UNIVERSE ));
    FORALLfacets {
      for (i = 0; i < dim; i++)
        facetv[i] = -facet->normal[i];
      facetv[i] = -facet->offset;
      //NOTE: Printing facet->toporient here for the [0,1]^3 unit cube test case (see TestBConvex.h) demonstrates that facet->toporient is meaningless in the qhull output.
      Constraint c = constraint_generator(convert_vector(facetv, facetq));
      ((C_Polyhedron*)m_poly)->add_constraint(c);
      num_facets++;
    }
    FORALLvertices {
      num_vertices++;
    }
#else
    m_poly = new_poly(dim);
    FORALLfacets {
      num_facets++;
    }
    FORALLvertices {
      for (i = 0; i < dim; i++)
        vertexv[i] = vertex->point[i];
      Generator g = point_generator(convert_vector(vertexv, vertexq));
      ((C_Polyhedron*)m_poly)->add_generator(g);
      num_vertices++;
    }
#endif
    //std::cout << "Qhull: " << num_facets << " facets and " << num_vertices << " vertices" << std::endl;
    //std::cout << "PPL: " << this->num_facets() << " facets and " << this->num_vertices() << " vertices" << std::endl;
    
    VW_ASSERT(!this->empty(), LogicErr() << "Convex hull is empty!");
    VW_ASSERT(!((const C_Polyhedron*)m_poly)->is_universe(), LogicErr() << "Convex hull is universe!");
    
    qhull_free();
#else // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
    return;
#endif // defined(VW_HAVE_PKG_QHULL) && VW_HAVE_PKG_QHULL==1
  }
  
  void BConvex::grow_( Vector<Promoted> const& point ) {
    Vector<mpq_class> pointq;
    ((C_Polyhedron*)m_poly)->add_generator(point_generator(convert_vector(point, pointq)));
  }

  void BConvex::grow( BConvex const& bconv ) {
    new_poly_if_needed(((const C_Polyhedron*)(bconv.poly()))->space_dimension(), m_poly);
    ((C_Polyhedron*)m_poly)->poly_hull_assign(*((const C_Polyhedron*)(bconv.poly())));
  }

  void BConvex::crop( BConvex const& bconv ) {
    if (!m_poly)
      return;
    ((C_Polyhedron*)m_poly)->intersection_assign(*((const C_Polyhedron*)(bconv.poly())));
  }

  void BConvex::expand( double offset ) {
    VW_ASSERT(!empty(), LogicErr() << "Cannot expand an empty polyhedron!");
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<mpq_class> center_( dim );
    poly_center((const C_Polyhedron*)m_poly, center_);
    void *poly = new_poly(dim);
    Vector<mpq_class> p( dim ), d;
    double norm;
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    Generator_System::const_iterator i;
    for (i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, p);
      d = p - center_;
      norm = norm_2_gmp(d);
      if (offset < 0.0 && -offset >= norm)
        p = center_;
      else
        p = center_ + (1.0 + offset / norm) * d;
      Generator g = point_generator(p);
      ((C_Polyhedron*)poly)->add_generator(g);
    }
    delete (C_Polyhedron*)m_poly;
    m_poly = poly;
  }

  bool BConvex::contains_( Vector<Promoted> const& point ) const {
    Vector<mpq_class> pointq;
    return (((const C_Polyhedron*)m_poly)->relation_with(point_generator(convert_vector(point, pointq))) == Poly_Gen_Relation::subsumes());
  }

  bool BConvex::contains( BConvex const& bconv ) const {
    if (!m_poly)
      return false;
    return ((const C_Polyhedron*)m_poly)->contains(*((const C_Polyhedron*)(bconv.poly())));
  }

  bool BConvex::intersects( BConvex const& bconv ) const {
    if (!m_poly)
      return false;
    return !((const C_Polyhedron*)m_poly)->is_disjoint_from(*((const C_Polyhedron*)(bconv.poly())));
  }

  double BConvex::size() const {
    if (!m_poly)
      return 0.0;
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    double result = 0;
    double d;
    Vector<mpq_class> v(dim);
    Vector<mpq_class> c;
    poly_center((const C_Polyhedron*)m_poly, c);
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, v);
      d = norm_2_sqr_gmp(v - c);
      result = std::max(result, d);
    }
    result = 2*std::sqrt(result);
    return result;
  }

  Vector<double> BConvex::center() const {
    VW_ASSERT(!empty(), LogicErr() << "Cannot find the center of an empty polyhedron!");
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<mpq_class> c(dim);
    Vector<double> cd(dim);
    poly_center((const C_Polyhedron*)m_poly, c);
    unconvert_vector(c, cd);
    return cd;
  }

  bool BConvex::equal( BConvex const& bconv ) const {
    if (empty())
      return bconv.empty();
    return *((const C_Polyhedron*)m_poly) == *((const C_Polyhedron*)(bconv.poly()));
  }

  bool BConvex::empty() const {
    if (!m_poly)
      return true;
    return ((const C_Polyhedron*)m_poly)->is_empty();
  }

  void BConvex::print( std::ostream& os ) const {
    if (!m_poly) {
      os << "false";
      return;
    }
    using ::Parma_Polyhedra_Library::IO_Operators::operator<<;
    os << *((const C_Polyhedron*)m_poly);
  }

  void BConvex::write( std::ostream& os, bool binary /* = false*/ ) const {
    if (!empty()) {
      unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
      std::vector<Vector<double> > points;
      Vector<mpq_class> v(dim);
      Vector<double> vd(dim);
      const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
      for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
        convert_point(*i, v);
        unconvert_vector(v, vd);
        points.push_back(vd);
      }
      write_point_list(os, points, binary);
    }
  }

  void BConvex::write( const char *fn, bool binary /* = false*/ ) const {
    if (binary) {
      std::ofstream of(fn, std::ofstream::out | std::ofstream::binary);
      write(of, binary);
      of.close();
    }
    else {
      std::ofstream of(fn);
      write(of, binary);
      of.close();
    }
  }

  void BConvex::write_vrml( std::ostream& os ) const {
    os << "#VRML V1.0 ascii\n\n";
    os << "# Created by the Intelligent Robotics Group,\n";
    os << "# NASA Ames Research Center\n";
    os << "# File generated by the NASA Ames Vision Workbench.\n\n";
    os << "Separator {\n";
  
    if (!empty()) {
      unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
      unsigned num_points;
      unsigned num_facets;
      Vector<double> vertexv( dim );
      std::vector<Vector<mpq_class> > points;
      std::vector<std::vector<unsigned> > facets;
      unsigned i, j, k;
      
      poly_facet_list((const C_Polyhedron*)m_poly, points, facets);
      num_points = points.size();
      num_facets = facets.size();
      
      os << "   Coordinate3 {\n";
      os << "      point [\n";
      for (i = 0; i < num_points; i++) {
        unconvert_vector(points[i], vertexv);
        os << "       ";
        for (k = 0; k < dim; k++)
          os << " " << vertexv[k];
        for (; k < 3; k++)
          os << " 0";
        os << ",\n";
      }
      os << "      ]\n";
      os << "   }\n";
      
      os << "   Material {\n";
      os << "      ambientColor 1 1 1\n";
      os << "      diffuseColor 1 1 1\n";
      os << "      specularColor 0.00 0.00 0.00\n";
      os << "      emissiveColor 0.00 0.00 0.00\n";
      os << "      shininess 0.00\n";
      os << "      transparency 0.00\n";
      os << "   }\n";
      os << "   IndexedFaceSet {\n";
      os << "      coordIndex [\n";
      os << "      ";
      for (i = 0; i < num_facets; i++) {
        os << "   " << facets[i][0] << ",";
        for (j = 1; j < facets[i].size(); j++)
          os << " " << facets[i][j] << ",";
        os << " " << facets[i][0] << ",";
        os << " -1,\n";
        os << "      ";
      }
      os << "]\n";
      os << "   }\n";
      
      os << "   Material {\n";
      os << "      ambientColor 0 0 1\n";
      os << "      diffuseColor 0 0 1\n";
      os << "      specularColor 0.00 0.00 0.00\n";
      os << "      emissiveColor 0.00 0.00 0.00\n";
      os << "      shininess 0.00\n";
      os << "      transparency 0.00\n";
      os << "   }\n";
      os << "   IndexedLineSet {\n";
      os << "      coordIndex [\n";
      os << "      ";
      for (i = 0; i < num_facets; i++) {
        os << "   " << facets[i][0] << ",";
        for (j = 1; j < facets[i].size(); j++)
          os << " " << facets[i][j] << ",";
        os << " " << facets[i][0] << ",";
        os << " -1,\n";
        os << "      ";
      }
      os << "]\n";
      os << "   }\n";
    }
    
    os << "}\n";
  }

  void BConvex::write_vrml( const char *fn ) const {
    std::ofstream of(fn);
    write_vrml(of);
  }
  
  void BConvex::write_oogl( std::ostream& os ) const {
    os << "OFF\n\n";
    os << "# Created by the Intelligent Robotics Group,\n";
    os << "# NASA Ames Research Center\n";
    os << "# File generated by the NASA Ames Vision Workbench.\n\n";
  
    if (empty()) {
      os << "0 0 0\n";
    }
    else {
      unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
      unsigned num_points;
      unsigned num_facets;
      unsigned num_edges;
      Vector<double> vertexv( dim );
      std::vector<Vector<mpq_class> > points;
      std::vector<std::vector<unsigned> > facets;
      unsigned i, j, k;
      
      poly_facet_list((const C_Polyhedron*)m_poly, points, facets);
      num_points = points.size();
      num_facets = facets.size();
      
      num_edges = 0;
      for (i = 0; i < num_facets; i++) {
        if (facets[i].size() > 1)
          num_edges += facets[i].size();
      }
      num_edges /= 2;
      
      os << num_points << " " << num_facets << " " << num_edges << "\n\n";
      
      for (i = 0; i < num_points; i++) {
        unconvert_vector(points[i], vertexv);
        for (k = 0; k < dim; k++)
          os << vertexv[k] << " ";
        for (; k < 3; k++)
          os << "0 ";
        os << "\n";
      }
      os << "\n";
      
      for (i = 0; i < num_facets; i++) {
        os << facets[i].size();
        for (j = 0; j < facets[i].size(); j++)
          os << " " << facets[i][j];
        os << " 1 1 1 1\n";
      }
    }
  }

  void BConvex::write_oogl( const char *fn ) const {
    std::ofstream of(fn);
    write_oogl(of);
  }
  
  unsigned BConvex::num_facets() const {
    if (!m_poly)
      return 0;
    return num_constraints(((const C_Polyhedron*)m_poly)->minimized_constraints());
  }

  unsigned BConvex::num_vertices() const {
    if (!m_poly)
      return 0;
    return num_generators(((const C_Polyhedron*)m_poly)->minimized_generators());
  }
  
  BBoxN BConvex::bounding_box() const {
    BBoxN bbox;
    if (!m_poly)
      return bbox;
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<mpq_class> v(dim);
    Vector<double> vd(dim);
    const Generator_System &gs = ((const C_Polyhedron*)m_poly)->minimized_generators();
    for (Generator_System::const_iterator i = gs.begin(); i != gs.end(); i++) {
      convert_point(*i, v);
      unconvert_vector(v, vd);
      bbox.grow(vd);
    }
    return bbox;
  }

  void BConvex::operator_mult_eq_( Promoted const& s ) {
    mpq_class q;
    convert_scalar(s, q);
    scale_poly(q, (C_Polyhedron*)m_poly);
  }

  void BConvex::operator_div_eq_( Promoted const& s ) {
    mpq_class q;
    convert_scalar(s, q);
    scale_poly(mpq_class(q.get_den(), q.get_num()), (C_Polyhedron*)m_poly);
  }

  void BConvex::operator_plus_eq_( Vector<Promoted> const& v ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<mpq_class> vq( dim );
    convert_vector(v, vq);
    offset_poly(vq, (C_Polyhedron*)m_poly);
  }

  void BConvex::operator_minus_eq_( Vector<Promoted> const& v ) {
    unsigned dim = ((const C_Polyhedron*)m_poly)->space_dimension();
    Vector<mpq_class> vq( dim );
    convert_vector(v, vq);
    offset_poly(-vq, (C_Polyhedron*)m_poly);
  }

  std::ostream& operator<<( std::ostream& os, BConvex const& bconv ) {
    bconv.print(os);
    return os;
  }

  void write_bconvex( std::string const& filename, BConvex const& bconv, bool binary /* = false*/ ) {
    bconv.write(filename.c_str(), binary);
  }

  void read_bconvex( std::string const& filename, BConvex& bconv, bool binary /* = false*/ ) {
    std::vector<Vector<double> > points;
    read_point_list(filename, points, binary);
    bconv = BConvex(points);
  }

}} // namespace vw::math
