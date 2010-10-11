// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ParticleSwarmOptimization.h
///
/// Non-linear Particle Swarm Optimization (PSO)
///
/// Based on: J. Kennedy and R. Eberhart, Particle Swarm Optimization, 1995
/// Digital Object Identifier: 10.1109/ICNN.1995.488968
/// http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=488968&isnumber=10434
///
/// Further reading: http://en.wikipedia.org/wiki/Particle_swarm_optimization
///
/// Running a PSO requires passing a functor object whose operator() accepts
/// an argument of type DomainT (e.g. a vw::Vector<double, n> and returns a
/// scalar value.
/// Also the upper and lower search domain boundaries - both of type DomainT
/// have to be supplied.
///
/// A number of arithmetic and other operations have to be defined on DomainT,
/// thus using vectors vw::Vector<double, n> or vw::Vector<double> is recommended.

#ifndef __VW_MATH_PARTICLE_SWARM_OPTIMIZATION_H__
#define __VW_MATH_PARTICLE_SWARM_OPTIMIZATION_H__

// Vision Workbench
#include <vw/Math/Functions.h>
#include <vw/Math/Vector.h>
#include <vw/Core/Log.h>

namespace vw {
namespace math {
  template <class FuncT, class DomainT>
  DomainT particle_swarm_optimization( FuncT const& func, DomainT const& min, DomainT const& max,
                                       bool verbose = false, int restarts = 1,
                                       unsigned int n_particles = 100,  unsigned int n_iter = 1000,
                                       double w = 0.9, double c1 = 2, double c2 = 2, double v_max = 4.0)
  {
    std::vector<DomainT> x;        x.resize(n_particles);         // particles
    std::vector<DomainT> x_hat;    x_hat.resize(n_particles);     // local maxima
    std::vector<double> x_hat_val; x_hat_val.resize(n_particles); // local maxima values
    std::vector<DomainT> v;        v.resize(n_particles);         // velocities

    std::srand(std::time(0)); // seed random number generator

    DomainT g_hat = min; // global minimum
    double g_hat_val = func(min);

    // external loop - run the whole PSO multiple times in case it gets in some local optimum
    // in function of the initial conditions and/or the development over time
    for (int start = 0; start < restarts; start++) {
      if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start+1 << "/" << restarts << std::endl;

      // initialize x random and v zero
      for (unsigned int i = 0; i < x.size(); i++) {
        x[i].set_size(min.size());
        x_hat[i].set_size(min.size());
        v[i].set_size(min.size());

        for (unsigned int j = 0; j < min.size(); j++) {
            x[i](j) = static_cast<double>(rand())/RAND_MAX*(max(j) - min(j)) + min(j);
            x_hat[i](j) = x[i](j);
            v[i](j) = 0;
        }

        x_hat_val[i] = func(x_hat[i]);
      }

      if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start+1 << "/" << restarts << " initialized particles" << std::endl;

      // search globally minimum particle
      for (unsigned int i = 0; i < x.size(); i++) {
        if (double d = func(x[i]) < g_hat_val) {
          g_hat_val = d;
          g_hat = x[i];
        }
      }

      if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start+1 << "/" << restarts << " found globally max particle with val " << g_hat_val << std::endl;

      // create two random number vectors with elements 0..1
      DomainT r1, r2;
      r1.set_size(min.size());
      r2.set_size(min.size());

      // iterate
      for (unsigned int iter = 0; iter < n_iter; iter++) {
        if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start << "/" << restarts << " iteration " << iter << "/" << n_iter << std::endl;

        // update all particles
        for (unsigned int i = 0; i < x.size(); i++) {
          // initialize random vectors
          for (unsigned int j = 0; j < r1.size(); j++) {
            r1(j) = static_cast<double>(rand())/RAND_MAX;
            r2(j) = static_cast<double>(rand())/RAND_MAX;
          }

          // particle position and velocity update
          x[i] = x[i] + v[i];
          v[i] = w*v[i] + c1*elem_prod(r1, x_hat[i] - x[i]) + c2*elem_prod(r2, g_hat - x[i]);

          // enforce maximum velocity
          double norm = vw::math::norm_2_sqr(v[i]);
          if (norm > v_max*v_max)
            v[i] /= vw::sqrt(norm);

          if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO x = " << x[i] << " with velocity v = " << v[i] << std::endl;

          // update local maximum
          double x_i_val = func(x[i]);
          if (x_i_val < x_hat_val[i]) {
            if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start+1 << "/" << restarts << " iteration " << iter << "/" << n_iter << " new local minimum " << x_i_val << std::endl;
            x_hat[i] = x[i];
            x_hat_val[i] = x_i_val;
          }

          // update global maximum
          if (x_i_val < g_hat_val) {
            if (verbose) vw::vw_out(vw::VerboseDebugMessage, "math") << "PSO run " << start+1 << "/" << restarts << " iteration " << iter << "/" << n_iter << " new global minimum " << x_i_val << std::endl;
            g_hat = x[i];
            g_hat_val = x_i_val;
          }
        }
      }
    }

    return g_hat;
  }
}
}

#endif /* __VW_MATH_PARTICLE_SWARM_OPTIMIZATION_H__ */

