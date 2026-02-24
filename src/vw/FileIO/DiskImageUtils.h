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
#include <vw/Image/ImageChannels.h>

namespace vw {

  /// Write a vector object to disk as an image file.  This function
  /// is particularly useful if you write the vector as an OpenEXR
  /// image file; this retains the maximal amount of information.
  template <class T>
  void write_vector(const std::string &filename, vw::Vector<T> &out_vector) {

    // Convert the vector to an image so that we can write it using
    // write_image().  There is probably a more efficient way to do
    // this, but this is the simplest way for now.
    vw::ImageView<T> out_image(out_vector.size(), 1, 1);

    for (size_t i = 0; i < out_vector.size(); i++)
      out_image(i, 0) = out_vector.impl()(i);

    write_image(filename, out_image);
  }

  /// Read a vector object from an image file on disk.  This function
  /// is particularly useful if the vector was saved as an OpenEXR
  /// image file; this retains the maximal amount of information.
  template <class T>
  void read_vector(vw::Vector<T>& in_vector, const std::string &filename) {

    // Treat the vector as an image so we can read from it using read_image().
    // There is probably a more efficient way to do this, but this is the
    // simplest way for now.
    vw::ImageView<T> buffer_image;
    read_image(buffer_image, filename);
    VW_ASSERT(buffer_image.planes() == 1,
              vw::IOErr() << "read_vector: Image file must be monochromatic"
              << " (1 plane/channel) to be read into a vector");
    vw::Vector<T> result(buffer_image.cols());

    for (int i = 0; i < buffer_image.cols(); i++)
      result.impl()(i) = buffer_image(i, 0);
    
    in_vector = result;
  }

  /// Find how many channels/bands are in a given image
  inline int get_num_channels(std::string filename){
    boost::shared_ptr<vw::DiskImageResource> src(vw::DiskImageResourcePtr(filename));
    return src->channels()*src->planes();
  }

  template <int m, int n, class T>
  struct SelectChannels: public vw::ReturnFixedType< vw::Vector<T, m>> {
    int k;
    SelectChannels(int k_in):k(k_in){}
    
    vw::Vector<T, m> operator() (vw::Vector<T, n> const& pt) const {
      return subvector(pt,k,m);
    }
  };

  /// Function that extracts the first m channels starting at channel k
  /// of an image with n channels. Must have m and n as templated arguments.
  template <int m, int n, class T>
  vw::UnaryPerPixelView<vw::DiskImageView< vw::Vector<T, n>>,
                        SelectChannels<m, n, T>>
  inline select_channels(vw::ImageViewBase<vw::DiskImageView< vw::Vector<T, n>>>
                          const& image, int k) {
    return vw::UnaryPerPixelView<vw::DiskImageView< vw::Vector<T, n>>,
      SelectChannels<m, n, T>>(image.impl(), SelectChannels<m, n, T>(k));
  }
  
  /// Function that extracts the first m channels starting at channel k
  /// of an image with n channels. Must have m and n as templated arguments.
  template <int m, int n, class T>
  vw::UnaryPerPixelView<vw::ImageViewRef< vw::Vector<T, n>>,
                        SelectChannels<m, n, T>>
  inline select_channels(vw::ImageViewBase<vw::ImageViewRef< vw::Vector<T, n>>>
                          const& image, int k) {
    return vw::UnaryPerPixelView<vw::ImageViewRef< vw::Vector<T, n>>,
      SelectChannels<m, n, T>>(image.impl(), SelectChannels<m, n, T>(k));
  }

  /// Read m channels from an image starting with channel k. Must
  /// have m as a templated argument. Channel starts at 0.
  /// Compile-intensive: results in 12 template instantiations.
  /// Prefer the wrappers in ImageChannelRead.h or PointCloudRead.h.
  template<int m, class T>
  vw::ImageViewRef<vw::Vector<T, m>> read_channels(std::string const& filename, int k) {
    
    int max_n = 12; // Maximum number of channels we can handle
    int n = get_num_channels(filename);
    
    VW_ASSERT(2 <= n, // This avoids a runtime error
            vw::ArgumentErr() << "Attempting to read an image with " << n
            << " channel(s) into an image with vector pixels.");
    VW_ASSERT(0 <= k,
               vw::ArgumentErr() << "Attempting to read channel " << k
               << " from an image.");
    VW_ASSERT(1 <= m,
               vw::ArgumentErr() << "Attempting to read " << m
               << " channel(s) from an image.");
    VW_ASSERT(k + m <= n,
               vw::ArgumentErr() << "Attempting to read " << k + m
               << " channels from an image with " << n << " channel(s).");
    VW_ASSERT(n <= max_n,
               vw::NoImplErr() << "Reading from images with more than "
               << max_n << " channels is not implemented.");

    // Turning off rescaling is important. Otherwise images with
    // uint16 channels cannot be read properly as floats.
    boost::shared_ptr<vw::DiskImageResource> rsrc(vw::DiskImageResourcePtr(filename));
    rsrc->set_rescale(false);
    
    vw::ImageViewRef<vw::Vector<T, m>> out_image;
    if      (n == 1) out_image = vw::select_channels<m, 1, T>
      (vw::DiskImageView<vw::Vector<T, 1>>(rsrc), k);
    else if (n == 2) out_image = vw::select_channels<m, 2, T>
      (vw::DiskImageView<vw::Vector<T, 2>>(rsrc), k);
    else if (n == 3) out_image = vw::select_channels<m, 3, T>
      (vw::DiskImageView<vw::Vector<T, 3>>(rsrc), k);
    else if (n == 4) out_image = vw::select_channels<m, 4, T>
      (vw::DiskImageView<vw::Vector<T, 4>>(rsrc), k);
    else if (n == 5) out_image = vw::select_channels<m, 5, T>
      (vw::DiskImageView<vw::Vector<T, 5>>(rsrc), k);
    else if (n == 6) out_image = vw::select_channels<m, 6, T>
      (vw::DiskImageView<vw::Vector<T, 6>>(rsrc), k);
    else if (n == 7) out_image = vw::select_channels<m, 7, T>
      (vw::DiskImageView<vw::Vector<T, 7>>(rsrc), k);
    else if (n == 8) out_image = vw::select_channels<m, 8, T>
      (vw::DiskImageView<vw::Vector<T, 8>>(rsrc), k);
    else if (n == 9) out_image = vw::select_channels<m, 9, T>
      (vw::DiskImageView<vw::Vector<T, 9>>(rsrc), k);
    else if (n == 10) out_image = vw::select_channels<m, 10, T>
      (vw::DiskImageView<vw::Vector<T, 10>>(rsrc), k);
    else if (n == 11) out_image = vw::select_channels<m, 11, T>
      (vw::DiskImageView<vw::Vector<T, 11>>(rsrc), k);
    else if (n == 12) out_image = vw::select_channels<m, 12, T>
      (vw::DiskImageView<vw::Vector<T, 12>>(rsrc), k);
    else
      vw::vw_throw(vw::NoImplErr() << "Reading from images with more than "
                   << max_n << " channels is not implemented.");

    return out_image;
  }

  // Read only one channel from an image. Channel starts at 0.
  // Compile-intensive via read_channels. Prefer ImageChannelRead.h.
  template<class T>
  vw::ImageViewRef<T> read_channel(std::string const& filename, int ch) {
    
    // Cannot read a Vector<T, 1> image, so have to handle this case separately.
    int n = vw::get_num_channels(filename);
    if (n > 1)
      return vw::select_channel(read_channels<1, T>(filename, ch), 0);
      
    return vw::DiskImageView<T>(filename);
  }
  
  // Load the first channel of an input image. Must load in its native
  // format then cast to double, in order to avoid auto-scaling of
  // pixels.
  vw::ImageViewRef<double> load_image_as_double(std::string const& image_file);

  // Returns the default nodata value for each data type (just the
  // minimum value).  Prefer to have
  // -std::numeric_limits<vw::float32>::max() as no-data, unless for
  // float64. Regardless of the inputs, the nodata value is read by vw
  // as a double anyway, and, for example, for uint8 data, setting 0
  // as no-data is not a good idea as it would conflict with valid
  // pixels.
  double get_default_nodata(int channel_type);
  
} // namespace vw

#endif
