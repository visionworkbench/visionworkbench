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

/// \file VectorIO.h
///
/// Functions for reading and writing vectors to files
///

//***
#ifndef __VW_MATH_VECTORIO_H__
#define __VW_MATH_VECTORIO_H__

#include <vw/Math/Vector.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResource.h>

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

    int i;
    for (i = 0; i < buffer_image.cols(); i++) {
	result.impl()(i) = buffer_image(i, 0);
    }
    in_vector = result;
  }

} // namespace vw

#endif // __VW_MATH_VECTORIO_H__


