// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// The helper source file is used to instantian different pixel
// versions of geoblend into different object files. It works by using
// the PIXEL_TYPE macro.

#include <vw/tools/geoblend.h>

using namespace vw;

#define INSTANTIATE_CUSTOM_BLEND(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_BLEND_(PIXELTYPEMACRO, CHANNELTYPE)

#define INSTANTIATE_CUSTOM_BLEND_(PIXELTYPEMACRO, CHANNELTYPE) INSTANTIATE_CUSTOM_BLEND__(PIXELTYPEMACRO, CHANNELTYPE, PIXELTYPEMACRO ## _ ## CHANNELTYPE)

#define INSTANTIATE_CUSTOM_BLEND__(PIXELTYPEMACRO, CHANNELTYPE, FUNCSUFFIX) \
  void do_blend_##FUNCSUFFIX(void) {                                    \
    do_blend<PIXELTYPEMACRO<CHANNELTYPE > >();                          \
  }

INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, uint8 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, int16 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, uint16 )
INSTANTIATE_CUSTOM_BLEND( PIXEL_TYPE, float32 )
