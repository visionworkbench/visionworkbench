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


/// \file Colormap.h
///
/// Algorithms for finding the colormap of an image

#ifndef __VW_IMAGE_COLORMAP_H__
#define __VW_IMAGE_COLORMAP_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>

#include <string>
#include <map>

namespace vw {
// These are specialized definitions for colormap functionality
  typedef Vector<uint8,3>                 Vector3u;
  typedef std::pair<std::string,Vector3u> lut_element;
  typedef std::vector<lut_element>        lut_type;

// Parse the colormap for given style. Note that if colormap_style
// is not one of the supported options, it can be a file name,
// from which the style is read. The output will be either in lut
// or lut_map, depending on the desired style. This will be sorted
// out later.
void parse_color_style(std::string const& colormap_style,
                       std::map<float, Vector3u> & lut_map);

// Populate a colormap from string pairs given in lut.
void populate_lut_map(lut_type const& lut,
                      std::map<float, Vector3u> & lut_map);

// Functor which colorizes every pixel
class Colormap: public ReturnFixedType<PixelMask<PixelRGB<uint8>>> {
  typedef std::map<float, Vector3u> map_type;
  map_type m_colormap;

public:
  Colormap(std::map<float, Vector3u> const& map) : m_colormap(map) {}

  PixelMask<PixelRGB<uint8>> operator()(PixelMask<PixelGray<float>> const& pix) const;
};

} // namespace vw

#endif // __VW_IMAGE_COLORMAP_H__
