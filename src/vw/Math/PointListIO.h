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

/// \file PointListIO.h
///
/// Provides functions to read/write point lists.
#ifndef __VW_MATH__POINT_LIST_IO_H__
#define __VW_MATH__POINT_LIST_IO_H__


#include <vector>
#include <iostream>
#include <fstream>

#include <vw/Math/Vector.h>
#include <vw/Core/Exception.h>


namespace vw {
namespace math {

  /// This function writes a list of points from any container that
  /// supports the size() and operator[] methods.  The container is
  /// usually a vw::Vector<>, but you could substitute other classes
  /// here as well.
  template <class ContainerT>
  void write_point_list(std::ostream &f, std::vector<ContainerT> const& pts) {
    VW_ASSERT( pts.size() > 0, LogicErr() << "No vectors to write!" );

    unsigned num_points = pts.size();
    unsigned dimensions = pts[0].size();

    f << dimensions << std::endl;
    f << num_points << std::endl;
    for (unsigned i = 0; i < num_points; i++) {
      for (unsigned j = 0; j < dimensions; j++)
        f << pts[i][j] << " ";
      f << std::endl;
    }
  }

  template <class ContainerT>
  void write_point_list(std::string const& filename, std::vector<ContainerT> const& pts) {
    VW_ASSERT( pts.size() > 0, LogicErr() << "No vectors to write!" );
    std::ofstream f(filename.c_str());
    VW_ASSERT( f, IOErr() << "Unable to open file for writing!" );
    write_point_list(f, pts);
    f.close();
  }

  /// This function reads a list of points into any container that
  /// supports the set_size() and operator[] methods.  The container is
  /// usually a vw::Vector<>, but you could substitute other classes
  /// here as well.
  template <class ContainerT>
  void read_point_list(std::istream &f, std::vector<ContainerT> &pts) {
    ContainerT v;
    double d;
    unsigned num_points;
    unsigned dimensions;

    f >> dimensions;
    VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
    f >> num_points;
    VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
    v.set_size(dimensions);
    for (unsigned i = 0; i < num_points; i++) {
      for (unsigned j = 0; j < dimensions; j++) {
        f >> d;
        VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
        v[j] = d;
      }
      pts.push_back(v);
    }
  }

  template <class ContainerT>
  void read_point_list(std::string const& filename, std::vector<ContainerT> &pts) {
    std::ifstream f(filename.c_str());
    VW_ASSERT( f, IOErr() << "Unable to open file for reading!" );
    read_point_list(f, pts);
    f.close();
  }

}} // namespace vw::math

#endif // __VW_MATH__POINT_LIST_IO_H__
