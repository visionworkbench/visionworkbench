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

/// \file NelderMead.h
/// 
/// Non-linear optimization using a simplex method.  
///
/// Based on Chapter 13, Section 13.1 of "The Nature of Mathematical
/// Modelling" by Neal Gershenfeld with some insprication from
/// Numerical Recipes.

#ifndef __VW_MATH_NELDER_MEAD_H__
#define __VW_MATH_NELDER_MEAD_H__

// Vision workbench
#include <vw/Core/Debugging.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

namespace vw {
namespace math {

  namespace optimization {
    enum NM_STATUS_CODES { eNelderMeadDidNotConverge = -1, 
                           eNelderMeadConvergedRelTolerance = 1 };
  }
  
  template <class FuncT, class DomainT>
  class Simplex {

    // Store the vertices in a pair: [ position, function value (cost) at that location ]
    typedef std::pair<DomainT, double> vertex_type;
    typedef typename std::list<vertex_type>::iterator vertex_iterator;

    FuncT m_func;
    std::list<vertex_type> m_vertices;

    // Insert a new vertex into the list while maintaining the
    // ordering such that f(v[1]) > f(v[2]) > ... > f(v[n])
    void insert_vertex(vertex_type vertex) {
      vertex_iterator iter = m_vertices.begin();

      while (vertex.second < (*iter).second)
        iter++;
      m_vertices.insert(iter, vertex);
    }
    
    // Search for the vertex in the list that had the highest function
    // evaluation value and return an iterator to that position.
    vertex_type& highest_vertex() { return m_vertices.front(); }
    vertex_type& lowest_vertex() { return m_vertices.back(); }

    vertex_type pop_highest_vertex() { 
      vertex_type vertex = m_vertices.front();
      m_vertices.pop_front();
      return vertex;
    }

    DomainT mean_vertex_location() { 
      vertex_iterator iter = m_vertices.begin();
      DomainT mean_location = (*iter).first;
      while ( ++iter != m_vertices.end())
        mean_location += (*iter).first;

      return mean_location / m_vertices.size();
    }


  public:
    template <class ScaleT>
    Simplex(FuncT const& func, DomainT const& seed, ScaleT const& scales) : m_func(func) {

      VW_ASSERT(scales.size() == seed.size(), 
                ArgumentErr() << "NelderMeadMinimizer: the number of scales does not match the dimensionality of the data in the seed vector.");

      // Create an N+1 vertex simplex where N is the dimensionality of
      // the domain.  
      m_vertices.clear();
      
      // The seed vector sets the center vertex of the simplex and the
      // other vertices are chosen by picking points in each of the
      // standard unit vector directions.  The distance in each
      // direction is determined by the scale variable.
      m_vertices.push_front( vertex_type(seed, m_func(seed)) );

      // Note -- this is probably not the most memory efficient way to
      // generate vectors from the identity matrix...
      vw::Matrix<double> identity = identity_matrix(seed.size());
      for (int i= 1; i < seed.size()+1; ++i) {
        DomainT vertex = seed + scales[i-1] * select_col(identity, i-1);
        insert_vertex( vertex_type(vertex, m_func(vertex)) ); 
      }
    }

    // Print out the current set of simplex vertices along with the
    // score at each vertex.
    void print_vertices() {
      for ( vertex_iterator iter = m_vertices.begin();
            iter != m_vertices.end(); ++iter) {
        std::cout << "\tSIMPLEX: " << (*iter).first << "[" << (*iter).second << "]   ";
      }
      std::cout << "\n";
    }

    DomainT location() { return lowest_vertex().first; }
    double value() { return lowest_vertex().second; }

    // This function causes the simplex to make one step down the
    // optimization surface.  Returns the size of the step that was
    // ultimately made.
    double update() {
      vertex_type highest_vtx = pop_highest_vertex();
      DomainT mean_vtx_loc = mean_vertex_location();
      vertex_type& lowest_vtx = lowest_vertex();

      // First, try to reflect the highest_vtx across the
      // mean_vtx. (reflect)
      DomainT new_loc = 2*mean_vtx_loc - highest_vtx.first;
      double new_val = m_func(new_loc); 

      // If this is now the lowest point, we should move another step
      // in this direction because this direction is a good
      // one. (reflect & extend)
      if (new_val < lowest_vtx.second) {
        DomainT reflect_extend_loc = 3*mean_vtx_loc - 2*highest_vtx.first;
        double reflect_extend_val = m_func(reflect_extend_loc); 

        if (reflect_extend_val < new_val) {
          new_loc = reflect_extend_loc;
          new_val = reflect_extend_val;
        }

      }

      // On the other hand, if the the reflected point is still the
      // highest point, we overshot the minimum so we will next try
      // reflecting and shrinking. (reflect & shrink)
      if (new_val > highest_vertex().second) {
        new_loc = 1.5*mean_vtx_loc - 0.5*highest_vtx.first;
        new_val = m_func(new_loc);
      }
      
      // If this is still worse, we try to shrinking the largest
      // vertex without reflecting it. (shrink)
      if (new_val > highest_vertex().second) {
        new_loc = 0.5 * (mean_vtx_loc + highest_vtx.first);
        new_val = m_func(new_loc);
      }

      // If this is _still_ worse, we assume that our step size is too
      // large so we shrink the entire simplex toward the best point
      // in the simplex (the best point stays in the same
      // place). (shrink all)
      if (new_val > highest_vertex().second) {
        new_loc = 0.5 * (highest_vtx.first + lowest_vtx.first);
        new_val = m_func(new_loc);

        vertex_iterator iter = m_vertices.begin();
        while(iter != m_vertices.end()) {
          // Omit the last position
          vertex_iterator next_iter = iter;
          if (++next_iter != m_vertices.end()) {
            (*iter).first = 0.5 * ((*iter).first + lowest_vtx.first);
            (*iter).second = m_func((*iter).first);
          }
          ++iter;
        }
      }

      // Finally, add in the moved vertex back into the simplex.
      insert_vertex( vertex_type(new_loc, new_val) );

      // ... and return the amount of improvement.
      return highest_vtx.second - new_val;
    }
  };

  template <class FuncT, class DomainT, class ScaleT>
  DomainT nelder_mead( FuncT const& func, DomainT const& seed, ScaleT const& scale, 
                       int &status, bool verbose = false, int restarts = 1, 
                       double tolerance = 1e-16, int max_iterations = 1000) {
    DomainT result = seed;
    status = optimization::eNelderMeadConvergedRelTolerance;

    if (verbose) {
      std::cout << "Nelder Mead Optimizer:\n";
      std::cout << "\tTol: " << tolerance << "   MaxIter: " << max_iterations << "   Restarts: " << restarts << "\n";
    }

    // Restart the simplex several times -- this prevents false
    // termination, which can happend from time to time if the simplex
    // gets stuck.
    int iterations;
    for (int i=0; i < restarts; ++i) {
      iterations = 0;
      Simplex<FuncT, DomainT> simplex(func, result, scale);
      
      // Perform simplex updates until tolerance in reached or
      // max_iterations is reached
      double delta = simplex.update();
      while (fabs(delta) > tolerance && iterations <= max_iterations) {
        delta = simplex.update();
        ++iterations;
        if (verbose && iterations % 100 == 0)
          std::cout << "\t" << iterations << ": " << simplex.location() << "[" << simplex.value() << "]\n";
      }

      if (verbose)
        std::cout << "\t" << iterations << ": " << simplex.location() << "[" << simplex.value() << "]\n";

      // Store the best result so far.
      result = simplex.location();
    }

    if (iterations >= max_iterations)
      status = optimization::eNelderMeadDidNotConverge;

    return result;
  }

  template <class FuncT, class DomainT>
    DomainT nelder_mead( FuncT const& func, DomainT const& seed, int &status, bool verbose = false,
                       int restarts = 1, double tolerance = 1e-16, int max_iterations = 1000) {
    vw::Vector<double> scale(seed.size());
    fill(scale,1.0);
    return nelder_mead(func, seed, scale, status, verbose, restarts, tolerance, max_iterations);
  }
  
}} // namespace vw::math

#endif // __VW_MATH_NELDER_MEAD_H__
