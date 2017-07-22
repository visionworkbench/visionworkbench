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


/// \file DiskImageUtils.h

#ifndef __VW_FILEIO_DISKIMAGEUTILS_H__
#define __VW_FILEIO_DISKIMAGEUTILS_H__

#include <vw/FileIO/DiskImageView.h>
#include <vw/Image/ImageViewRef.h>

namespace vw {

  /// Write a vector object to disk as an image file.  This function
  /// is particularly useful if you write the vector as an OpenEXR
  /// image file; this retains the maximal amount of information.
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

  /// Read a vector object from an image file on disk.  This function
  /// is particularly useful if the vector was saved as an OpenEXR
  /// image file; this retains the maximal amount of information.
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

  /// Find how many channels/bands are in a given image
  inline int get_num_channels(std::string filename){
    boost::shared_ptr<vw::DiskImageResource> src(vw::DiskImageResourcePtr(filename));
    return src->channels()*src->planes();
  }

  template <int m, int n, class T>
  struct SelectChannels : public vw::ReturnFixedType< vw::Vector<T, m> > {
    int k;
    SelectChannels(int k_in):k(k_in){}
    
    vw::Vector<T, m> operator() (vw::Vector<T, n> const& pt) const {
      return subvector(pt,k,m);
    }
  };

  /// Function that extracts the first m channels starting at channel k
  /// of an image with n channels. Must have m and n as templated arguments.
  template <int m, int n, class T>
  vw::UnaryPerPixelView<vw::DiskImageView< vw::Vector<T, n> >,
                        SelectChannels<m, n, T> >
  inline select_channels( vw::ImageViewBase<vw::DiskImageView< vw::Vector<T, n> > >
                          const& image, int k ) {
    return vw::UnaryPerPixelView<vw::DiskImageView< vw::Vector<T, n> >,
      SelectChannels<m, n, T> >( image.impl(), SelectChannels<m, n, T>(k) );
  }
  
  /// Function that extracts the first m channels starting at channel k
  /// of an image with n channels. Must have m and n as templated arguments.
  template <int m, int n, class T>
  vw::UnaryPerPixelView<vw::ImageViewRef< vw::Vector<T, n> >,
                        SelectChannels<m, n, T> >
  inline select_channels( vw::ImageViewBase<vw::ImageViewRef< vw::Vector<T, n> > >
                          const& image, int k ) {
    return vw::UnaryPerPixelView<vw::ImageViewRef< vw::Vector<T, n> >,
      SelectChannels<m, n, T> >( image.impl(), SelectChannels<m, n, T>(k) );
  }

  /// Read m channels from an image starting with channel k. Must
  /// have m as a templated argument.
  template<int m, class T>
  vw::ImageViewRef< vw::Vector<T, m> > read_channels(std::string const& filename,
                                                     int k){
    
    int max_n = 6;
    int n = get_num_channels(filename);

    VW_ASSERT( 0 <= k,
               vw::ArgumentErr() << "Attempting to read channel " << k
               << " from an image.");
    VW_ASSERT( 1 <= m,
               vw::ArgumentErr() << "Attempting to read " << m
               << " channel(s) from an image.");
    VW_ASSERT( k + m <= n,
               vw::ArgumentErr() << "Attempting to read " << k + m
               << " channels from an image with " << n << " channel(s).");
    VW_ASSERT( n <= max_n,
               vw::NoImplErr() << "Reading from images with more than "
               << max_n << " channels is not implemented.");

    // Turning off rescaling is important. Otherwise images with
    // uint16 channels cannot be read properly as floats.
    boost::shared_ptr<vw::DiskImageResource> rsrc( vw::DiskImageResourcePtr(filename));
    rsrc->set_rescale(false);
    
    vw::ImageViewRef< vw::Vector<T, m> > out_image;
    if      (n == 1) out_image = select_channels<m, 1, T>
      (vw::DiskImageView< vw::Vector<T, 1> >(rsrc), k);
    else if (n == 2) out_image = select_channels<m, 2, T>
      (vw::DiskImageView< vw::Vector<T, 2> >(rsrc), k);
    else if (n == 3) out_image = select_channels<m, 3, T>
      (vw::DiskImageView< vw::Vector<T, 3> >(rsrc), k);
    else if (n == 4) out_image = select_channels<m, 4, T>
      (vw::DiskImageView< vw::Vector<T, 4> >(rsrc), k);
    else if (n == 5) out_image = select_channels<m, 5, T>
      (vw::DiskImageView< vw::Vector<T, 5> >(rsrc), k);
    else if (n == 6) out_image = select_channels<m, 6, T>
      (vw::DiskImageView< vw::Vector<T, 6> >(rsrc), k);

    return out_image;
  }
  
 } // namespace vw

#endif
