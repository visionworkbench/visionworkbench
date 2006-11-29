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

/// \file GlobalTonemap.h
/// 
/// Tonemapping operators that operate on the entire image.  That is,
/// the tonemapped value of each pixel can be computed without knowing
/// the values of pixels in its immediate neighborhood.
///
/// The tonemapping operators currently defined in this file include:
///
/// - Drago Logarithmic Tonemap Operator
/// 

#ifndef __VW_HDR_GLOBALTONEMAP_H__
#define __VW_HDR_GLOBALTONEMAP_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

namespace vw { 
namespace hdr {

  ImageView<PixelRGB<double> > drago_tone_map(ImageView<PixelRGB<double> > hdr_image, 
                                              double bias = 0.85, 
                                              double exposure_factor = 1.0, 
                                              double L_dmax = 100);

}} // namespace vw::HDR


#endif  // __GlobalToneMap_H__
