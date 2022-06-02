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

/// \file MatrixIO.cc
///
/// Functions for reading and writing matrices to files
///
#include <vw/FileIO/MatrixIO.h>

namespace vw {

  // Write a double precision matrix to disk stored in plain text.
  void write_matrix_as_txt(std::string const& filename, vw::Matrix<double> const& matrix) {
  std::ofstream ofs(filename.c_str());
    ofs.precision(17);
    
    for (int row = 0; row < matrix.rows(); row++) {
      for (int col = 0; col < matrix.cols(); col++) {
        ofs << matrix(row, col);
        if (col + 1 < matrix.cols()) 
          ofs << " ";
        else
          ofs << "\n";
      }
    }

    return;
  }
  
  void read_matrix_as_txt(std::string const& filename, vw::Matrix<double> & matrix) {

    // First read the values from each row
    std::ifstream ifs(filename.c_str());
    std::vector<std::vector<double>> mat;
    std::string line;
    while (getline(ifs, line)){
      std::istringstream iss(line.c_str());
      std::vector<double> vec;
      double val = -1.0;
      while (iss >> val)
        vec.push_back(val);

      mat.push_back(vec);
    }

    // See if all rows have the same number of values
    for (size_t it = 1; it < mat.size(); it++) {
      if (mat[0].size() != mat[it].size()) {
        vw::vw_throw( ArgumentErr() << "Not all rows have the same number of values in: "
                      << filename << "\n" );
      }
    }

    matrix.set_size(mat.size(), mat[0].size());
    for (int row = 0; row < matrix.rows(); row++) {
      for (int col = 0; col < matrix.cols(); col++) {
        matrix(row, col) = mat[row][col];
      }
    }
  }
  
} // namespace vw
