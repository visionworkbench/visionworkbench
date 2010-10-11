// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file LocalToneMap.h
///
/// Tonemapping operators that compute the tonemapping in one small
/// image region at a time.
///
/// This file implements the following tone mapping operators.
///
///
#ifndef __VW_HDR_LOCALTONEMAP_H__
#define __VW_HDR_LOCALTONEMAP_H__

#include <vw/Image/ImageView.h>

namespace vw {
namespace hdr {

  ImageView<PixelRGB<double> > ashikhmin_tone_map(ImageView<PixelRGB<double> > hdr_image,
                                                  double threshold = 0.5);

}} // namespace vw::HDR

#endif  // __VW_HDR_LOCALTONEMAP_H__
