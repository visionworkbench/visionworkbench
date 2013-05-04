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


/// \file Palette.h
///
/// Simple support for color palettes.  A lot more should be here
/// that currentlly isn't, including a palette data type, the ability
/// to convert from true color to palette image, and tools to pick an
/// optimal color palette for an image.
///
#ifndef __VW_IMAGE_PALETTE_H__
#define __VW_IMAGE_PALETTE_H__

#include <vw/Core/Functors.h>
#include <vw/Image/ImageView.h>

namespace vw {

  // A palette filter.  Loads in a palette stored as an
  // image file.  Converts floating-point pixels in the
  // range [0,1] to the apprpriate values by indexing into
  // the (first row of the) image.  Yes, it's a hack.
  template <class PixelT>
  class PaletteFilter : public ReturnFixedType<PixelT>
  {
    ImageView<PixelT> palette;
    int entries;
  public:
    PaletteFilter( std::string const& palette_file ) {
      read_image( palette, palette_file );
      entries = palette.cols();
    }

    PixelT operator()( float val ) const {
      int index = lround( val * entries );
      if( index >= entries ) index = entries-1;
      else if( index < 0 ) index = 0;
      return palette( index, 0 );
    }
  };

} // namespace vw

#endif // __VW_IMAGE_PALETTE_H__
