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

/// \file PixelTypeInfo.cc
/// 
/// Functions that return information about types by ID.
///
#include <vw/Image/PixelTypeInfo.h>

vw::int32 vw::channel_size( vw::ChannelTypeEnum type ) {
  switch( type ) {
  case VW_CHANNEL_BOOL:
    return sizeof(bool);
  case VW_CHANNEL_CHAR:
  case VW_CHANNEL_INT8:
  case VW_CHANNEL_UINT8:
  case VW_CHANNEL_GENERIC_1_BYTE:
    return 1;
  case VW_CHANNEL_INT16:
  case VW_CHANNEL_UINT16:
  case VW_CHANNEL_GENERIC_2_BYTE:
    return 2;
  case VW_CHANNEL_INT32:
  case VW_CHANNEL_UINT32:
  case VW_CHANNEL_FLOAT32:
  case VW_CHANNEL_GENERIC_4_BYTE:
    return 4;
  case VW_CHANNEL_INT64:
  case VW_CHANNEL_UINT64:
  case VW_CHANNEL_FLOAT64:
  case VW_CHANNEL_GENERIC_8_BYTE:
    return 8;
  default:
    vw_throw( ArgumentErr() << "Unrecognized or unsupported channel type (" << type << ")." );
    return 0; // never reached
  }
}

const char *vw::channel_type_name( vw::ChannelTypeEnum format ) {
  switch( format ) {
  case VW_CHANNEL_BOOL: return "BOOL";
  case VW_CHANNEL_CHAR: return "CHAR";
  case VW_CHANNEL_INT8: return "INT8";
  case VW_CHANNEL_UINT8: return "UINT8";
  case VW_CHANNEL_INT16: return "INT16";
  case VW_CHANNEL_UINT16: return "UINT16";
  case VW_CHANNEL_INT32: return "INT32";
  case VW_CHANNEL_UINT32: return "UINT32";
  case VW_CHANNEL_FLOAT32: return "FLOAT32";
  case VW_CHANNEL_INT64: return "INT64";
  case VW_CHANNEL_UINT64: return "UINT64";
  case VW_CHANNEL_FLOAT64: return "FLOAT64";
  default: return "UNKNOWN";
  }
}

vw::int32 vw::num_channels( vw::PixelFormatEnum format ) {
  switch( format ) {
  case VW_PIXEL_SCALAR:
  case VW_PIXEL_GRAY:
  case VW_PIXEL_GENERIC_1_CHANNEL:
    return 1;
  case VW_PIXEL_GRAYA:
  case VW_PIXEL_GENERIC_2_CHANNEL:
    return 2;
  case VW_PIXEL_RGB:
  case VW_PIXEL_HSV:
  case VW_PIXEL_XYZ:
  case VW_PIXEL_GENERIC_3_CHANNEL:
    return 3;
  case VW_PIXEL_RGBA:
  case VW_PIXEL_GENERIC_4_CHANNEL:
    return 4;
  default:
    vw_throw( ArgumentErr() << "Unrecognized or unsupported pixel format (" << format << ")." );
    return 0; // never reached
  }
}

const char *vw::pixel_format_name( vw::PixelFormatEnum format ) {
  switch( format ) {
  case VW_PIXEL_SCALAR: return "SCALAR";
  case VW_PIXEL_GRAY: return "GRAY";
  case VW_PIXEL_GRAYA: return "GRAYA";
  case VW_PIXEL_RGB: return "RGB";
  case VW_PIXEL_RGBA: return "RGBA";
  case VW_PIXEL_HSV: return "HSV";
  case VW_PIXEL_XYZ: return "XYZ";
  default: return "UNKNOWN";
  }
}
