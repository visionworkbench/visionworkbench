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

/// \file MatrixIO.h
///
/// Functions for reading and writing matrices to files
///

//***
#ifndef __VW_FILEIO_MATRIXIO_H__
#define __VW_FILEIO_MATRIXIO_H__

#include <vw/Math/Matrix.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageIO.h>

//Matrix<double> alignMatrix;
//write_matrix(out_prefix + "-align.exr", align_matrix);

namespace vw {

  /// Write a matrix object to disk as an image file.  This function
  /// is particularly useful if you write the matrix as an OpenEXR
  /// image file; this retains the maximal amount of information.
  template <class T>
  void write_matrix(const std::string &filename, vw::Matrix<T> &out_matrix) {

    // Convert the matrix to an image so that we can write it using
    // write_image().  There is probably a more efficient way to do
    // this, but this is the simplest way for now.
    vw::ImageView<T> out_image(out_matrix.cols(), out_matrix.rows(), 1);

    unsigned int i, j;
    for (i = 0; i < out_matrix.cols(); i++) {
      for (j = 0; j < out_matrix.rows(); j++) {
        out_image(i, j) = out_matrix.impl()(j, i);
      }
    }
    write_image(filename, out_image);
  }


  /// Read a matrix object from an image file on disk.  This function
  /// is particularly useful if the matrix was saved as an OpenEXR
  /// image file; this retains the maximal amount of information.
  template <class T>
  void read_matrix(vw::Matrix<T>& in_matrix, const std::string &filename) {

    // Treat the matrix as an image so we can read from it using read_image().
    // There is probably a more efficient way to do this, but this is the
    // simplest way for now.
    vw::ImageView<T> buffer_image;
    read_image(buffer_image, filename);
    VW_ASSERT(buffer_image.planes() == 1,
              vw::IOErr() << "read_matrix: Image file must be monochromatic"
              << " (1 plane/channel) to be read into a matrix");
    vw::Matrix<T> result(buffer_image.rows(), buffer_image.cols());

    int i, j;
    for (i = 0; i < buffer_image.cols(); i++) {
      for (j = 0; j < buffer_image.rows(); j++) {
        result.impl()(j, i) = buffer_image(i, j);
      }
    }
    in_matrix = result;
  }

  // Write a double precision matrix to disk stored in plain text.
  void write_matrix_as_txt(std::string const& filename, vw::Matrix<double> const& matrix);

  // Read a double precision matrix to disk stored in plain text.
  void read_matrix_as_txt(std::string const& filename, vw::Matrix<double> & matrix);
  
} // namespace vw

#endif // __VW_FILEIO_MATRIXIO_H__


