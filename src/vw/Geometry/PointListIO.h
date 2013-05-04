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


/// \file PointListIO.h
///
/// Provides functions to read/write point lists.
#ifndef __VW_MATH__POINT_LIST_IO_H__
#define __VW_MATH__POINT_LIST_IO_H__


#include <vector>
#include <iostream>
#include <fstream>

#include <vw/Core/Exception.h>


namespace vw {
namespace math {

  /// This function writes a list of points from any container that
  /// supports the size() and operator[] methods.  The container is
  /// usually a vw::Vector<>, but you could substitute other classes
  /// here as well.
  template <class ContainerT>
  void write_point_list(std::ostream &f, std::vector<ContainerT> const& pts, bool binary = false) {
    VW_ASSERT( pts.size() > 0, LogicErr() << "No vectors to write!" );

    unsigned num_points = pts.size();
    unsigned dimensions = pts[0].size();

    if (binary) {
      f.write((char*)&dimensions, sizeof(dimensions));
      f.write((char*)&num_points, sizeof(num_points));
      for (unsigned i = 0; i < num_points; i++) {
        for (unsigned j = 0; j < dimensions; j++)
          f.write((char*)&pts[i][j], sizeof(pts[i][j]));
      }
    }
    else {
      f << dimensions << std::endl;
      f << num_points << std::endl;
      for (unsigned i = 0; i < num_points; i++) {
        for (unsigned j = 0; j < dimensions; j++)
          f << pts[i][j] << " ";
        f << std::endl;
      }
    }
  }

  template <class ContainerT>
  void write_point_list(std::string const& filename, std::vector<ContainerT> const& pts, bool binary = false) {
    VW_ASSERT( pts.size() > 0, LogicErr() << "No vectors to write!" );

    if (binary) {
      std::ofstream f(filename.c_str(), std::ofstream::out | std::ofstream::binary);
      VW_ASSERT( f, IOErr() << "Unable to open file for writing!" );
      write_point_list(f, pts, binary);
      f.close();
    }
    else {
      std::ofstream f(filename.c_str());
      VW_ASSERT( f, IOErr() << "Unable to open file for writing!" );
      write_point_list(f, pts, binary);
      f.close();
    }
  }

  /// This function reads a list of points into any container that
  /// supports the set_size() and operator[] methods.  The container is
  /// usually a vw::Vector<>, but you could substitute other classes
  /// here as well.
  template <class ContainerT>
  void read_point_list(std::istream &f, std::vector<ContainerT> &pts, bool binary = false) {
    ContainerT v;
    unsigned num_points;
    unsigned dimensions;

    if (binary) {
      f.read((char*)&dimensions, sizeof(dimensions));
      VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
      f.read((char*)&num_points, sizeof(num_points));
      VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
      v.set_size(dimensions);
      for (unsigned i = 0; i < num_points; i++) {
        for (unsigned j = 0; j < dimensions; j++) {
          f.read((char*)&v[j], sizeof(v[j]));
          VW_ASSERT( !f.fail(), IOErr() << "Invalid point list file format!" );
        }
        pts.push_back(v);
      }
    }
    else {
      double d;
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
  }

  template <class ContainerT>
  void read_point_list(std::string const& filename, std::vector<ContainerT> &pts, bool binary = false) {
    if (binary) {
      std::ifstream f(filename.c_str(), std::ifstream::in | std::ifstream::binary);
      VW_ASSERT( f, IOErr() << "Unable to open file for reading!" );
      read_point_list(f, pts, binary);
      f.close();
    }
    else {
      std::ifstream f(filename.c_str());
      VW_ASSERT( f, IOErr() << "Unable to open file for reading!" );
      read_point_list(f, pts, binary);
      f.close();
    }
  }

}} // namespace vw::math

#endif // __VW_MATH__POINT_LIST_IO_H__
