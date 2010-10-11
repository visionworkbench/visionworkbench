// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file MatrixIO.h
///
/// Functions for reading and writing matrices to files
///

//***
#ifndef __VW_INTERESTPOINT_MATRIXIO_H__
#define __VW_INTERESTPOINT_MATRIXIO_H__

#include <vw/Math/Matrix.h>
#include <vw/Image/ImageView.h>

//Matrix<double> alignMatrix;
//write_matrix(out_prefix + "-align.exr", align_matrix);

namespace vw {

  // Write a matrix object to disk as an image file.  This function
  // is particularly useful if you write the matrix as an OpenEXR
  // image file; this retains the maximal amount of information.
  template <class T>
  void write_matrix( const std::string &filename, vw::Matrix<T> &out_matrix ) {

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


  // Read a matrix object from an image file on disk.  This function
  // is particularly useful if the matrix was saved as an OpenEXR
  // image file; this retains the maximal amount of information.
  template <class T>
  void read_matrix( vw::Matrix<T>& in_matrix, const std::string &filename ) {

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

} // namespace vw

#endif // __VW_INTERESTPOINT_MATRIXIO_H__


