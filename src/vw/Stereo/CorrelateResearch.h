// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_CORRELATE_RESEARCH_H__
#define __VW_STEREO_CORRELATE_RESEARCH_H__

#include <vw/Stereo/Correlate.h>

namespace vw {
namespace stereo {

  template <class ChannelT>
  void subpixel_correlation_affine_2d(ImageView<PixelMask<Vector2f> > &disparity_map,
                                      ImageView<ChannelT> const& left_image,
                                      ImageView<ChannelT> const& right_image,
                                      int kern_width, int kern_height,
                                      bool do_horizontal_subpixel = true,
                                      bool do_vertical_subpixel = true,
                                      bool verbose = false);

  template <class ChannelT>
  void subpixel_correlation_affine_2d_bayesian(ImageView<PixelMask<Vector2f> > &disparity_map,
                                               ImageView<ChannelT> const& left_image,
                                               ImageView<ChannelT> const& right_image,
                                               int kern_width, int kern_height,
                                               bool do_horizontal_subpixel = true,
                                               bool do_vertical_subpixel = true,
                                               bool verbose = false);


}}

#endif//__VW_STEREO_CORRELATE_RESEARCH_H__
