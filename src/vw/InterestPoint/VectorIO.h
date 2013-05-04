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


/// \file VectorIO.h
///
/// Functions for reading and writing vectors to files
///

//***
#ifndef __VW_MATH_VECTORIO_H__
#define __VW_MATH_VECTORIO_H__

#include <vw/Math/Vector.h>
#include <vw/Image/ImageView.h>

//Vector<double> alignVector;
//write_vector(out_prefix + "-align.exr", align_vector);

namespace vw {

  // Write a vector object to disk as an image file.  This function
  // is particularly useful if you write the vector as an OpenEXR
  // image file; this retains the maximal amount of information.
  template <class T>
  void write_vector( const std::string &filename, vw::Vector<T> &out_vector ) {

    // Convert the vector to an image so that we can write it using
    // write_image().  There is probably a more efficient way to do
    // this, but this is the simplest way for now.
    vw::ImageView<T> out_image(out_vector.size(), 1, 1);

    unsigned int i;
    for (i = 0; i < out_vector.size(); i++) {
      out_image(i, 0) = out_vector.impl()(i);
    }
    write_image(filename, out_image);
  }


  // Read a vector object from an image file on disk.  This function
  // is particularly useful if the vector was saved as an OpenEXR
  // image file; this retains the maximal amount of information.
  template <class T>
  void read_vector( vw::Vector<T>& in_vector, const std::string &filename ) {

    // Treat the vector as an image so we can read from it using read_image().
    // There is probably a more efficient way to do this, but this is the
    // simplest way for now.
    vw::ImageView<T> buffer_image;
    read_image(buffer_image, filename);
    VW_ASSERT(buffer_image.planes() == 1,
              vw::IOErr() << "read_vector: Image file must be monochromatic"
              << " (1 plane/channel) to be read into a vector");
    vw::Vector<T> result(buffer_image.cols());

    for (int i = 0; i < buffer_image.cols(); i++) {
      result.impl()(i) = buffer_image(i, 0);
    }
    in_vector = result;
  }

} // namespace vw

#endif // __VW_MATH_VECTORIO_H__


