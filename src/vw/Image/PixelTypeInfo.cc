// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PixelTypeInfo.cc
///
/// Functions that return information about types by ID.
///
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Core/Exception.h>

#include <boost/algorithm/string.hpp>

bool vw::simple_conversion(vw::ChannelTypeEnum a, vw::ChannelTypeEnum b) {
  if (a == b)
    return true;
  if (channel_size_nothrow(a) != channel_size_nothrow(b))
    return false;

  switch (a) {
    case VW_CHANNEL_GENERIC_1_BYTE:
    case VW_CHANNEL_GENERIC_2_BYTE:
    case VW_CHANNEL_GENERIC_4_BYTE:
    case VW_CHANNEL_GENERIC_8_BYTE: return true;
    default: /* noop */ break;
  }

  switch (b) {
    case VW_CHANNEL_GENERIC_1_BYTE:
    case VW_CHANNEL_GENERIC_2_BYTE:
    case VW_CHANNEL_GENERIC_4_BYTE:
    case VW_CHANNEL_GENERIC_8_BYTE: return true;
    default: /* noop */ break;
  }

  // we can do something about CHAR <-> UINT8/INT8, but it's rarely used, and we
  // have to figure out if char is signed or unsigned.

  return false;
}

bool vw::simple_conversion(vw::PixelFormatEnum a, vw::PixelFormatEnum b) {
  if (a == b)
    return true;
  if (num_channels_nothrow(a) != num_channels_nothrow(b))
    return false;

  switch (a) {
    case VW_PIXEL_SCALAR:
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GENERIC_1_CHANNEL:
    case VW_PIXEL_GENERIC_2_CHANNEL:
    case VW_PIXEL_GENERIC_3_CHANNEL:
    case VW_PIXEL_GENERIC_4_CHANNEL:
    case VW_PIXEL_GENERIC_5_CHANNEL:
    case VW_PIXEL_GENERIC_6_CHANNEL: return true;
    default: /* noop */ break;
  }

  switch (b) {
    case VW_PIXEL_SCALAR:
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GENERIC_1_CHANNEL:
    case VW_PIXEL_GENERIC_2_CHANNEL:
    case VW_PIXEL_GENERIC_3_CHANNEL:
    case VW_PIXEL_GENERIC_4_CHANNEL:
    case VW_PIXEL_GENERIC_5_CHANNEL:
    case VW_PIXEL_GENERIC_6_CHANNEL: return true;
    default: /* noop */ break;
  }

  return false;
}

vw::uint32 vw::channel_size_nothrow( vw::ChannelTypeEnum type ) {
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
    case VW_CHANNEL_FLOAT16:
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
      return 0;
  }
}

vw::uint32 vw::channel_size( vw::ChannelTypeEnum type ) {
  vw::uint32 size = channel_size_nothrow(type);
  VW_ASSERT(size > 0, ArgumentErr() << "Unrecognized or unsupported channel type (" << type << ")." );
  return size;
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
  case VW_CHANNEL_FLOAT16: return "FLOAT16";
  case VW_CHANNEL_FLOAT32: return "FLOAT32";
  case VW_CHANNEL_INT64: return "INT64";
  case VW_CHANNEL_UINT64: return "UINT64";
  case VW_CHANNEL_FLOAT64: return "FLOAT64";
  case VW_CHANNEL_GENERIC_1_BYTE: return "GENERIC_1_BYTE";
  case VW_CHANNEL_GENERIC_2_BYTE: return "GENERIC_2_BYTE";
  case VW_CHANNEL_GENERIC_4_BYTE: return "GENERIC_4_BYTE";
  case VW_CHANNEL_GENERIC_8_BYTE: return "GENERIC_8_BYTE";
  default: return "UNKNOWN";
  }
}

vw::uint32 vw::num_channels_nothrow( vw::PixelFormatEnum format ) {
  switch( format ) {
  case VW_PIXEL_SCALAR:
  case VW_PIXEL_GRAY:
  case VW_PIXEL_GENERIC_1_CHANNEL:
    return 1;
  case VW_PIXEL_GRAYA:
  case VW_PIXEL_SCALAR_MASKED:
  case VW_PIXEL_GRAY_MASKED:
  case VW_PIXEL_GENERIC_2_CHANNEL:
    return 2;
  case VW_PIXEL_RGB:
  case VW_PIXEL_HSV:
  case VW_PIXEL_XYZ:
  case VW_PIXEL_LUV:
  case VW_PIXEL_LAB:
  case VW_PIXEL_GRAYA_MASKED:
  case VW_PIXEL_GENERIC_3_CHANNEL:
    return 3;
  case VW_PIXEL_RGBA:
  case VW_PIXEL_RGB_MASKED:
  case VW_PIXEL_HSV_MASKED:
  case VW_PIXEL_XYZ_MASKED:
  case VW_PIXEL_LUV_MASKED:
  case VW_PIXEL_LAB_MASKED:
  case VW_PIXEL_GENERIC_4_CHANNEL:
    return 4;
  case VW_PIXEL_RGBA_MASKED:
  case VW_PIXEL_GENERIC_5_CHANNEL:
    return 5;
  case VW_PIXEL_GENERIC_6_CHANNEL:
    return 6;
  default:
    return 0;
  }
}

vw::uint32 vw::num_channels( vw::PixelFormatEnum format ) {
  vw::uint32 num = num_channels_nothrow(format);
  VW_ASSERT(num > 0, ArgumentErr() << "Unrecognized or unsupported pixel format (" << format << ")." );
  return num;
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
  case VW_PIXEL_LUV: return "LUV";
  case VW_PIXEL_LAB: return "LAB";
  case VW_PIXEL_UNKNOWN_MASKED: return "UNKNOWN_MASKED";
  case VW_PIXEL_SCALAR_MASKED: return "SCALAR_MASKED";
  case VW_PIXEL_GRAY_MASKED: return "GRAY_MASKED";
  case VW_PIXEL_GRAYA_MASKED: return "GRAYA_MASKED";
  case VW_PIXEL_RGB_MASKED: return "RGB_MASKED";
  case VW_PIXEL_RGBA_MASKED: return "RGBA_MASKED";
  case VW_PIXEL_HSV_MASKED: return "HSV_MASKED";
  case VW_PIXEL_XYZ_MASKED: return "XYZ_MASKED";
  case VW_PIXEL_LUV_MASKED: return "LUV_MASKED";
  case VW_PIXEL_LAB_MASKED: return "LAB_MASKED";
  case VW_PIXEL_GENERIC_1_CHANNEL: return "VW_PIXEL_GENERIC_1_CHANNEL";
  case VW_PIXEL_GENERIC_2_CHANNEL: return "VW_PIXEL_GENERIC_2_CHANNEL";
  case VW_PIXEL_GENERIC_3_CHANNEL: return "VW_PIXEL_GENERIC_3_CHANNEL";
  case VW_PIXEL_GENERIC_4_CHANNEL: return "VW_PIXEL_GENERIC_4_CHANNEL";
  case VW_PIXEL_GENERIC_5_CHANNEL: return "VW_PIXEL_GENERIC_5_CHANNEL";
  case VW_PIXEL_GENERIC_6_CHANNEL: return "VW_PIXEL_GENERIC_6_CHANNEL";
  default: return "UNKNOWN";
  }
}

vw::ChannelTypeEnum vw::channel_name_to_enum( const std::string& name ) {
  std::string uname(name);
  boost::to_upper(uname);

  if (uname.length() < 4 || uname.length() > 15)
    return VW_CHANNEL_UNKNOWN;

  switch(uname[0]) {
    case 'B':
      if (uname == "BOOL") return VW_CHANNEL_BOOL;
      break;
    case 'C':
      if (uname == "CHAR") return VW_CHANNEL_CHAR;
      break;
    case 'D':
      if (uname == "DOUBLE") return VW_CHANNEL_FLOAT64;
      break;
    case 'F':
      if (uname == "FLOAT16") return VW_CHANNEL_FLOAT16;
      if (uname == "FLOAT64") return VW_CHANNEL_FLOAT64;
      if (uname == "FLOAT32" || uname == "FLOAT") return VW_CHANNEL_FLOAT32;
      break;
    case 'G':
      if (uname == "GENERIC_1_BYTE") return VW_CHANNEL_GENERIC_1_BYTE;
      if (uname == "GENERIC_2_BYTE") return VW_CHANNEL_GENERIC_2_BYTE;
      if (uname == "GENERIC_4_BYTE") return VW_CHANNEL_GENERIC_4_BYTE;
      if (uname == "GENERIC_8_BYTE") return VW_CHANNEL_GENERIC_8_BYTE;
      break;
    case 'I':
      if (uname == "INT8")  return VW_CHANNEL_INT8;
      if (uname == "INT16") return VW_CHANNEL_INT16;
      if (uname == "INT32" || uname == "INT") return VW_CHANNEL_INT32;
      if (uname == "INT64") return VW_CHANNEL_INT64;
      break;
    case 'U':
      if (uname == "UINT8")  return VW_CHANNEL_UINT8;
      if (uname == "UINT16") return VW_CHANNEL_UINT16;
      if (uname == "UINT32" || uname == "UINT") return VW_CHANNEL_UINT32;
      if (uname == "UINT64") return VW_CHANNEL_UINT64;
      break;
    default:
      break;
  }
  return VW_CHANNEL_UNKNOWN;
}
