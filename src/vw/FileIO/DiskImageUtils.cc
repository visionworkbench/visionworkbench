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


/// \file DiskImageUtils.cc

#include <vw/FileIO/DiskImageUtils.h>

namespace vw {

// Load the first channel of an input image. Must load in its native
// format then cast to double, in order to avoid auto-scaling of
// pixels.
ImageViewRef<double> load_image_as_double(std::string const& image_file) {
  
  // Determining the format of the input images
  boost::shared_ptr<vw::DiskImageResource> rsrc(vw::DiskImageResourcePtr(image_file));
  ChannelTypeEnum input_data_type = rsrc->channel_type();

  ImageViewRef<double> image;
  
  switch (input_data_type) {
  case VW_CHANNEL_INT8:
    image = vw::pixel_cast<double>(DiskImageView<vw::int8>(image_file));    break;
  case VW_CHANNEL_UINT8:
    image = vw::pixel_cast<double>(DiskImageView<vw::uint8>(image_file));   break;
  case VW_CHANNEL_INT16:
    image = vw::pixel_cast<double>(DiskImageView<vw::int16>(image_file));   break;
  case VW_CHANNEL_UINT16:
    image = vw::pixel_cast<double>(DiskImageView<vw::uint16>(image_file));  break;
  case VW_CHANNEL_INT32:
    image = vw::pixel_cast<double>(DiskImageView<vw::int32>(image_file));   break;
  case VW_CHANNEL_UINT32:
    image = vw::pixel_cast<double>(DiskImageView<vw::uint32>(image_file));  break;
  case VW_CHANNEL_FLOAT32:
    image = vw::pixel_cast<double>(DiskImageView<vw::float32>(image_file)); break;
  case VW_CHANNEL_FLOAT64:
    image = vw::pixel_cast<double>(DiskImageView<vw::float64>(image_file)); break;
  default:
    vw_throw(ArgumentErr() << "Input image format " << input_data_type << " is not supported.\n");
  };

  return image;
}

// See .h for the documentation
double get_default_nodata(int channel_type) {
  switch (channel_type) {
  case    VW_CHANNEL_INT8   : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_UINT8  : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_INT16  : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_UINT16 : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_INT32  : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_UINT32 : return -std::numeric_limits<vw::float32>::max();
  case    VW_CHANNEL_FLOAT32: return -std::numeric_limits<vw::float32>::max();
  default:                    return -std::numeric_limits<vw::float64>::max(); // double
  };
}

} // namespace vw

