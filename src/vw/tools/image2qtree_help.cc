// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// The helper source file is used to instantian different pixel
// versions of image2qtree into different object files. It works by using
// the PIXEL_TYPE macro.

#include <vw/tools/image2qtree.h>

#define INSTANTIATE_CUSTOM_MOSAIC(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_MOSAIC_(PIXELTYPEMACRO, CHANNELTYPE)

#define INSTANTIATE_CUSTOM_MOSAIC_(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_MOSAIC__(PIXELTYPEMACRO, CHANNELTYPE, PIXELTYPEMACRO ## _ ## CHANNELTYPE)

#define INSTANTIATE_CUSTOM_MOSAIC__(PIXELTYPEMACRO, CHANNELTYPE, FUNCSUFFIX) \
  void do_mosaic_##FUNCSUFFIX(const Options& opt,                       \
                              const vw::ProgressCallback *progress) {   \
    do_mosaic<vw::PIXELTYPEMACRO<vw::CHANNELTYPE > >(opt, progress);    \
  }

INSTANTIATE_CUSTOM_MOSAIC( PIXEL_TYPE, CHANNEL_TYPE )
