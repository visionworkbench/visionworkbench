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

/// \file ImageOctaveHistory.h
/// 
/// Class for storing all of the intermediate images processed while
/// iterating through an ImageOctave. This can be useful for generating
/// descriptors after interest point detection has been completed.
/// 
#ifndef __IMAGE_OCTAVE_HISTORY_H__
#define __IMAGE_OCTAVE_HISTORY_H__

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/ImageOctave.h>

namespace vw {
namespace ip {

template <class ImageT>
class ImageOctaveHistory : std::vector<std::vector<ImageT> > {
 private:
  int num_scales;

 public:
  ImageOctaveHistory() : std::vector<std::vector<ImageT> >(0), num_scales(0) {}

  inline int octaves() const { return this->size(); }

  inline int scales() const {
    return num_scales;
  }

  inline void add_octave(const std::vector<ImageT>& octave) {
    this->push_back(octave);
    num_scales = octave.size() - 2;
  }

  const ImageT& image_at_scale(float scale) const {
    int octave = (int)(log(scale) / M_LN2);
    if (octave == octaves()) octave = octaves() - 1;
    int plane = ImageOctave<float>::scale_to_plane_index(1 << octave, num_scales, scale);
    return (*this)[octave][plane];
  }
};

}} // namespace vw::ip 

#endif // __IMAGE_OCTAVE_HISTORY_H__
