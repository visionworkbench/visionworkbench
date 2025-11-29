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


/// \file ImageResource.cc
///
/// Defines a run-type-typed image buffer.
///
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/BBox.h>
#include <vw/Image/ImageResource.h>

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#include <vector>
#endif
#include <map>
#include <cmath>

#include <boost/integer_traits.hpp>
#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/smart_ptr/shared_array.hpp>

using namespace vw;

// -----------------------------------------------------------------
// Section for handling type conversion


/// Declare function type: Convert a src value to a dest value
typedef void (*channel_convert_func)(void* src, void* dest);

/// Apply a cast from the src to the dest types
template <class SrcT, class DestT>
void channel_convert_cast( SrcT* src, DestT* dest ) {
  *dest = DestT(*src);
}

// Simple integer type conversion functions
void channel_convert_uint16_to_uint8( uint16* src, uint8* dest ) {
  *dest = uint8( *src / (65535/255) );
}
void channel_convert_uint8_to_uint16( uint8* src, uint16* dest ) {
  *dest = uint16( *src ) * (65535/255);
}

/// Convert any integer into a float in the -1 to +1 range.
/// First convert from int to double, do the multiplication, then cast to float.
template <class SrcT, class DestT>
void channel_convert_int_to_float(SrcT* src, DestT* dest) {
  *dest = DestT(double(*src) * (double(1.0)/double(boost::integer_traits<SrcT>::const_max)));
}

/// Convert a float in the range -1 to +1 to an integer type.
template <class SrcT, class DestT>
void channel_convert_float_to_int( SrcT* src, DestT* dest ) {
  if (*src > SrcT(1.0))
    *dest = boost::integer_traits<DestT>::const_max;
  else if (*src < SrcT(0.0))
    *dest = DestT(0);
  else
    *dest = DestT(double(*src) * double(boost::integer_traits<DestT>::const_max));
}

/// Pointers to two maps:  <type pair> -> conversion function
/// - One is for rescaling conversions, the other for non-rescaling.
std::map<std::pair<ChannelTypeEnum,ChannelTypeEnum>,channel_convert_func>
  *channel_convert_map = 0,
  *channel_convert_rescale_map = 0;

/// Helper class for adding entries to the two conversion maps mentioned above.
class ChannelConvertMapEntry {
private:
  /// Instantiate the two channel conversion maps mentioned above if they do not exist.
  void initialize() {
    if( !channel_convert_map )
      channel_convert_map = new std::map<std::pair<ChannelTypeEnum,ChannelTypeEnum>,channel_convert_func>();
    if( !channel_convert_rescale_map )
      channel_convert_rescale_map = new std::map<std::pair<ChannelTypeEnum,ChannelTypeEnum>,channel_convert_func>();
  }
public:
  /// Load a non-scaling conversion function into the two maps (the same function goes in both).
  template <class SrcT, class DstT>
  ChannelConvertMapEntry( void (*func)(SrcT*,DstT*) ) {
    initialize();
    ChannelTypeEnum src = ChannelTypeID<SrcT>::value;
    ChannelTypeEnum dst = ChannelTypeID<DstT>::value;
    channel_convert_map        ->operator[]( std::make_pair(src,dst) ) = (channel_convert_func)func;
    channel_convert_rescale_map->operator[]( std::make_pair(src,dst) ) = (channel_convert_func)func;
  }
  /// Load both a non-scaling and rescaling function into the two maps.
  template <class SrcT, class DstT>
  ChannelConvertMapEntry( void (*func)(SrcT*,DstT*), void (*rescale_func)(SrcT*,DstT*) ) {
    initialize();
    ChannelTypeEnum src = ChannelTypeID<SrcT>::value;
    ChannelTypeEnum dst = ChannelTypeID<DstT>::value;
    channel_convert_map        ->operator[]( std::make_pair(src,dst) ) = (channel_convert_func)func;
    channel_convert_rescale_map->operator[]( std::make_pair(src,dst) ) = (channel_convert_func)rescale_func;
  }
};

// Load up a bunch of generic conversion functions!
// - Might be better to use a default handler for the most common case.
ChannelConvertMapEntry _conv_i8i8  ( &channel_convert_cast<int8,int8>   );
ChannelConvertMapEntry _conv_i8u8  ( &channel_convert_cast<int8,uint8>  );
ChannelConvertMapEntry _conv_i8i16 ( &channel_convert_cast<int8,int16>  );
ChannelConvertMapEntry _conv_i8u16 ( &channel_convert_cast<int8,uint16> );
ChannelConvertMapEntry _conv_i8i32 ( &channel_convert_cast<int8,int32>  );
ChannelConvertMapEntry _conv_i8u32 ( &channel_convert_cast<int8,uint32> );
ChannelConvertMapEntry _conv_i8i64 ( &channel_convert_cast<int8,int64>  );
ChannelConvertMapEntry _conv_i8u64 ( &channel_convert_cast<int8,uint64> );
ChannelConvertMapEntry _conv_i8f32 ( &channel_convert_cast<int8,float >, &channel_convert_int_to_float<int8,float>  );
ChannelConvertMapEntry _conv_i8f64 ( &channel_convert_cast<int8,double>, &channel_convert_int_to_float<int8,double> );
ChannelConvertMapEntry _conv_u8i8  ( &channel_convert_cast<uint8,int8>   );
ChannelConvertMapEntry _conv_u8u8  ( &channel_convert_cast<uint8,uint8>  );
ChannelConvertMapEntry _conv_u8i16 ( &channel_convert_cast<uint8,int16>  );
ChannelConvertMapEntry _conv_u8u16 ( &channel_convert_cast<uint8,uint16>, &channel_convert_uint8_to_uint16 );
ChannelConvertMapEntry _conv_u8i32 ( &channel_convert_cast<uint8,int32>  );
ChannelConvertMapEntry _conv_u8u32 ( &channel_convert_cast<uint8,uint32> );
ChannelConvertMapEntry _conv_u8i64 ( &channel_convert_cast<uint8,int64>  );
ChannelConvertMapEntry _conv_u8u64 ( &channel_convert_cast<uint8,uint64> );
ChannelConvertMapEntry _conv_u8f32 ( &channel_convert_cast<uint8,float >, &channel_convert_int_to_float<uint8,float>  );
ChannelConvertMapEntry _conv_u8f64 ( &channel_convert_cast<uint8,double>, &channel_convert_int_to_float<uint8,double> );
ChannelConvertMapEntry _conv_i16i8 ( &channel_convert_cast<int16,int8>   );
ChannelConvertMapEntry _conv_i16u8 ( &channel_convert_cast<int16,uint8>  );
ChannelConvertMapEntry _conv_i16i16( &channel_convert_cast<int16,int16>  );
ChannelConvertMapEntry _conv_i16u16( &channel_convert_cast<int16,uint16> );
ChannelConvertMapEntry _conv_i16i32( &channel_convert_cast<int16,int32>  );
ChannelConvertMapEntry _conv_i16u32( &channel_convert_cast<int16,uint32> );
ChannelConvertMapEntry _conv_i16i64( &channel_convert_cast<int16,int64>  );
ChannelConvertMapEntry _conv_i16u64( &channel_convert_cast<int16,uint64> );
ChannelConvertMapEntry _conv_i16f32( &channel_convert_cast<int16,float >, &channel_convert_int_to_float<int16,float>  );
ChannelConvertMapEntry _conv_i16f64( &channel_convert_cast<int16,double>, &channel_convert_int_to_float<int16,double> );
ChannelConvertMapEntry _conv_u16i8 ( &channel_convert_cast<uint16,int8>   );
ChannelConvertMapEntry _conv_u16u8 ( &channel_convert_cast<uint16,uint8>, &channel_convert_uint16_to_uint8 );
ChannelConvertMapEntry _conv_u16i16( &channel_convert_cast<uint16,int16>  );
ChannelConvertMapEntry _conv_u16u16( &channel_convert_cast<uint16,uint16> );
ChannelConvertMapEntry _conv_u16i32( &channel_convert_cast<uint16,int32>  );
ChannelConvertMapEntry _conv_u16u32( &channel_convert_cast<uint16,uint32> );
ChannelConvertMapEntry _conv_u16i64( &channel_convert_cast<uint16,int64>  );
ChannelConvertMapEntry _conv_u16u64( &channel_convert_cast<uint16,uint64> );
ChannelConvertMapEntry _conv_u16f32( &channel_convert_cast<uint16,float >, &channel_convert_int_to_float<uint16,float>  );
ChannelConvertMapEntry _conv_u16f64( &channel_convert_cast<uint16,double>, &channel_convert_int_to_float<uint16,double> );
ChannelConvertMapEntry _conv_i32i8 ( &channel_convert_cast<int32,int8>   );
ChannelConvertMapEntry _conv_i32u8 ( &channel_convert_cast<int32,uint8>  );
ChannelConvertMapEntry _conv_i32i16( &channel_convert_cast<int32,int16>  );
ChannelConvertMapEntry _conv_i32u16( &channel_convert_cast<int32,uint16> );
ChannelConvertMapEntry _conv_i32i32( &channel_convert_cast<int32,int32>  );
ChannelConvertMapEntry _conv_i32u32( &channel_convert_cast<int32,uint32> );
ChannelConvertMapEntry _conv_i32i64( &channel_convert_cast<int32,int64>  );
ChannelConvertMapEntry _conv_i32u64( &channel_convert_cast<int32,uint64> );
ChannelConvertMapEntry _conv_i32f32( &channel_convert_cast<int32,float >, &channel_convert_int_to_float<int32,float>  );
ChannelConvertMapEntry _conv_i32f64( &channel_convert_cast<int32,double>, &channel_convert_int_to_float<int32,double> );
ChannelConvertMapEntry _conv_u32i8 ( &channel_convert_cast<uint32,int8>   );
ChannelConvertMapEntry _conv_u32u8 ( &channel_convert_cast<uint32,uint8>  );
ChannelConvertMapEntry _conv_u32i16( &channel_convert_cast<uint32,int16>  );
ChannelConvertMapEntry _conv_u32u16( &channel_convert_cast<uint32,uint16> );
ChannelConvertMapEntry _conv_u32i32( &channel_convert_cast<uint32,int32>  );
ChannelConvertMapEntry _conv_u32u32( &channel_convert_cast<uint32,uint32> );
ChannelConvertMapEntry _conv_u32i64( &channel_convert_cast<uint32,int64>  );
ChannelConvertMapEntry _conv_u32u64( &channel_convert_cast<uint32,uint64> );
ChannelConvertMapEntry _conv_u32f32( &channel_convert_cast<uint32,float >, &channel_convert_int_to_float<uint32,float>  );
ChannelConvertMapEntry _conv_u32f64( &channel_convert_cast<uint32,double>, &channel_convert_int_to_float<uint32,double> );
ChannelConvertMapEntry _conv_i64i8 ( &channel_convert_cast<int64,int8>   );
ChannelConvertMapEntry _conv_i64u8 ( &channel_convert_cast<int64,uint8>  );
ChannelConvertMapEntry _conv_i64i16( &channel_convert_cast<int64,int16>  );
ChannelConvertMapEntry _conv_i64u16( &channel_convert_cast<int64,uint16> );
ChannelConvertMapEntry _conv_i64i32( &channel_convert_cast<int64,int32>  );
ChannelConvertMapEntry _conv_i64u32( &channel_convert_cast<int64,uint32> );
ChannelConvertMapEntry _conv_i64i64( &channel_convert_cast<int64,int64>  );
ChannelConvertMapEntry _conv_i64u64( &channel_convert_cast<int64,uint64> );
ChannelConvertMapEntry _conv_i64f32( &channel_convert_cast<int64,float >, &channel_convert_int_to_float<int64,float>  );
ChannelConvertMapEntry _conv_i64f64( &channel_convert_cast<int64,double>, &channel_convert_int_to_float<int64,double> );
ChannelConvertMapEntry _conv_u64i8 ( &channel_convert_cast<uint64,int8>   );
ChannelConvertMapEntry _conv_u64u8 ( &channel_convert_cast<uint64,uint8>  );
ChannelConvertMapEntry _conv_u64i16( &channel_convert_cast<uint64,int16>  );
ChannelConvertMapEntry _conv_u64u16( &channel_convert_cast<uint64,uint16> );
ChannelConvertMapEntry _conv_u64i32( &channel_convert_cast<uint64,int32>  );
ChannelConvertMapEntry _conv_u64u32( &channel_convert_cast<uint64,uint32> );
ChannelConvertMapEntry _conv_u64i64( &channel_convert_cast<uint64,int64>  );
ChannelConvertMapEntry _conv_u64u64( &channel_convert_cast<uint64,uint64> );
ChannelConvertMapEntry _conv_u64f32( &channel_convert_cast<uint64,float >, &channel_convert_int_to_float<uint64,float > );
ChannelConvertMapEntry _conv_u64f64( &channel_convert_cast<uint64,double>, &channel_convert_int_to_float<uint64,double> );
ChannelConvertMapEntry _conv_f32i8 ( &channel_convert_cast<float, int8 >, &channel_convert_float_to_int<float, int8>  );
ChannelConvertMapEntry _conv_f32u8 ( &channel_convert_cast<float,uint8 >, &channel_convert_float_to_int<float,uint8>  );
ChannelConvertMapEntry _conv_f32i16( &channel_convert_cast<float, int16>, &channel_convert_float_to_int<float, int16> );
ChannelConvertMapEntry _conv_f32u16( &channel_convert_cast<float,uint16>, &channel_convert_float_to_int<float,uint16> );
ChannelConvertMapEntry _conv_f32i32( &channel_convert_cast<float, int32>, &channel_convert_float_to_int<float, int32> );
ChannelConvertMapEntry _conv_f32u32( &channel_convert_cast<float,uint32>, &channel_convert_float_to_int<float,uint32> );
ChannelConvertMapEntry _conv_f32i64( &channel_convert_cast<float, int64>, &channel_convert_float_to_int<float, int64> );
ChannelConvertMapEntry _conv_f32u64( &channel_convert_cast<float,uint64>, &channel_convert_float_to_int<float,uint64> );
ChannelConvertMapEntry _conv_f32f32( &channel_convert_cast<float,float > );
ChannelConvertMapEntry _conv_f32f64( &channel_convert_cast<float,double> );
ChannelConvertMapEntry _conv_f64i8 ( &channel_convert_cast<double, int8 >, &channel_convert_float_to_int<double, int8>  );
ChannelConvertMapEntry _conv_f64u8 ( &channel_convert_cast<double,uint8 >, &channel_convert_float_to_int<double,uint8>  );
ChannelConvertMapEntry _conv_f64i16( &channel_convert_cast<double, int16>, &channel_convert_float_to_int<double, int16> );
ChannelConvertMapEntry _conv_f64u16( &channel_convert_cast<double,uint16>, &channel_convert_float_to_int<double,uint16> );
ChannelConvertMapEntry _conv_f64i32( &channel_convert_cast<double, int32>, &channel_convert_float_to_int<double, int32> );
ChannelConvertMapEntry _conv_f64u32( &channel_convert_cast<double,uint32>, &channel_convert_float_to_int<double,uint32> );
ChannelConvertMapEntry _conv_f64i64( &channel_convert_cast<double, int64>, &channel_convert_float_to_int<double, int64> );
ChannelConvertMapEntry _conv_f64u64( &channel_convert_cast<double,uint64>, &channel_convert_float_to_int<double,uint64> );
ChannelConvertMapEntry _conv_f64f32( &channel_convert_cast<double,float > );
ChannelConvertMapEntry _conv_f64f64( &channel_convert_cast<double,double> );

//------------------------------------------------------------------------------------
// Section for assigning max value

// Each of these sections of code can be replaced with simpler functions
//  incorporating a switch statement: ChannelTypeEnum -> function call

// Function type: Assigns a channel the maximum value
typedef void (*channel_set_max_func)(void* dest);

// Implementation for any integer
template <class DestT>
void channel_set_max_int( DestT* dest ) {
  *dest = boost::integer_traits<DestT>::const_max;
}
// Implementation for any float type
template <class DestT>
void channel_set_max_float( DestT* dest ) {
  *dest = DestT(1.0);
}

// The map, class, and entry setting mirror the previous section.

std::map<ChannelTypeEnum,channel_set_max_func> *channel_set_max_map = 0;

class ChannelSetMaxMapEntry {
public:
  template <class DstT>
  ChannelSetMaxMapEntry( void (*func)(DstT*) ) {
    if( !channel_set_max_map )
      channel_set_max_map = new std::map<ChannelTypeEnum,channel_set_max_func>();
    ChannelTypeEnum dst = ChannelTypeID<DstT>::value;
    channel_set_max_map->operator[]( dst ) = (channel_set_max_func)func;
  }
};

ChannelSetMaxMapEntry _setmax_i8 ( &channel_set_max_int<int8  > );
ChannelSetMaxMapEntry _setmax_u8 ( &channel_set_max_int<uint8 > );
ChannelSetMaxMapEntry _setmax_i16( &channel_set_max_int<int16 > );
ChannelSetMaxMapEntry _setmax_u16( &channel_set_max_int<uint16> );
ChannelSetMaxMapEntry _setmax_i32( &channel_set_max_int<int32 > );
ChannelSetMaxMapEntry _setmax_u32( &channel_set_max_int<uint32> );
ChannelSetMaxMapEntry _setmax_i64( &channel_set_max_int<int64 > );
ChannelSetMaxMapEntry _setmax_u64( &channel_set_max_int<uint64> );
ChannelSetMaxMapEntry _setmax_f32( &channel_set_max_float<float > );
ChannelSetMaxMapEntry _setmax_f64( &channel_set_max_float<double> );

//------------------------------------------------------------------------------------
// Section for averaging channels

// Channel Average:
//   Reduces a number of channels into one by averaging
typedef void (*channel_average_func)(void* src, void* dest, int32 len);

template <class T>
void channel_average( T* src, T* dest, int32 len ) {
  typename AccumulatorType<T>::type accum = typename AccumulatorType<T>::type();
  for( int32 i=0; i<len; ++i ) accum += src[i];
  *dest = accum / len;
}

// The map, class, and entry setting mirror the previous section.

std::map<ChannelTypeEnum,channel_average_func> *channel_average_map = 0;

class ChannelAverageMapEntry {
public:
  template <class T>
  ChannelAverageMapEntry( void (*func)(T*,T*,int32) ) {
    if( !channel_average_map )
      channel_average_map = new std::map<ChannelTypeEnum,channel_average_func>();
    ChannelTypeEnum ctid = ChannelTypeID<T>::value;
    channel_average_map->operator[]( ctid ) = (channel_average_func)func;
  }
};

ChannelAverageMapEntry _average_i8 ( &channel_average<int8> );
ChannelAverageMapEntry _average_u8 ( &channel_average<uint8> );
ChannelAverageMapEntry _average_i16( &channel_average<int16> );
ChannelAverageMapEntry _average_u16( &channel_average<uint16> );
ChannelAverageMapEntry _average_i32( &channel_average<int32> );
ChannelAverageMapEntry _average_u32( &channel_average<uint32> );
ChannelAverageMapEntry _average_i64( &channel_average<int64> );
ChannelAverageMapEntry _average_u64( &channel_average<uint64> );
ChannelAverageMapEntry _average_f32( &channel_average<float> );
ChannelAverageMapEntry _average_f64( &channel_average<double> );


//------------------------------------------------------------------------------------
// Premultiply section

// Channel Premultiply:
//   Applies the Alpha Channel to the rest of the channels:
typedef void (*channel_premultiply_func)(void* src, void* dst, int32 len);

template <class T>
void channel_premultiply_int( T* src, T* dst, int32 len ) {
  double scale = src[len-1] / (double)(boost::integer_traits<T>::const_max);
  for( int32 i=0; i<len-1; ++i ) dst[i] = T( round(src[i] * scale) );
  dst[len-1] = src[len-1];
}

template <class T>
void channel_premultiply_float( T* src, T* dst, int32 len ) {
  double scale = (double)(src[len-1]);
  for( int32 i=0; i<len-1; ++i ) dst[i] = T( src[i] * scale );
  dst[len-1] = src[len-1];
}

// The map, class, and entry setting mirror the previous section.

std::map<ChannelTypeEnum,channel_premultiply_func> *channel_premultiply_map = 0;

class ChannelPremultiplyMapEntry {
public:
  template <class T>
  ChannelPremultiplyMapEntry( void (*func)(T*,T*,int32) ) {
    if( !channel_premultiply_map )
      channel_premultiply_map = new std::map<ChannelTypeEnum,channel_premultiply_func>();
    ChannelTypeEnum ctid = ChannelTypeID<T>::value;
    channel_premultiply_map->operator[]( ctid ) = (channel_premultiply_func)func;
  }
};

ChannelPremultiplyMapEntry _premultiply_i8 ( &channel_premultiply_int<int8> );
ChannelPremultiplyMapEntry _premultiply_u8 ( &channel_premultiply_int<uint8> );
ChannelPremultiplyMapEntry _premultiply_i16( &channel_premultiply_int<int16> );
ChannelPremultiplyMapEntry _premultiply_u16( &channel_premultiply_int<uint16> );
ChannelPremultiplyMapEntry _premultiply_i32( &channel_premultiply_int<int32> );
ChannelPremultiplyMapEntry _premultiply_u32( &channel_premultiply_int<uint32> );
ChannelPremultiplyMapEntry _premultiply_i64( &channel_premultiply_int<int64> );
ChannelPremultiplyMapEntry _premultiply_u64( &channel_premultiply_int<uint64> );
ChannelPremultiplyMapEntry _premultiply_f32( &channel_premultiply_float<float> );
ChannelPremultiplyMapEntry _premultiply_f64( &channel_premultiply_float<double> );

//-----------------------------------------------------------------
// Unpremultiply section

/// Channel Unpremultiply:
///   Removes the premultiply of alpha to other channels:
typedef void (*channel_unpremultiply_func)(void* src, void* dst, int32 len);

template <class T>
void channel_unpremultiply_int( T* src, T* dst, int32 len ) {
  double scale = src[len-1] / (double)(boost::integer_traits<T>::const_max);
  for( int32 i=0; i<len-1; ++i ) dst[i] = T( round(src[i] / scale) );
  dst[len-1] = src[len-1];
}

template <class T>
void channel_unpremultiply_float( T* src, T* dst, int32 len ) {
  double scale = (double)(src[len-1]);
  for( int32 i=0; i<len-1; ++i ) dst[i] = T( src[i] / scale );
  dst[len-1] = src[len-1];
}

// The map, class, and entry setting mirror the previous section.

std::map<ChannelTypeEnum,channel_unpremultiply_func> *channel_unpremultiply_map = 0;

class ChannelUnpremultiplyMapEntry {
public:
  template <class T>
  ChannelUnpremultiplyMapEntry( void (*func)(T*,T*,int32) ) {
    if( !channel_unpremultiply_map )
      channel_unpremultiply_map = new std::map<ChannelTypeEnum,channel_unpremultiply_func>();
    ChannelTypeEnum ctid = ChannelTypeID<T>::value;
    channel_unpremultiply_map->operator[]( ctid ) = (channel_unpremultiply_func)func;
  }
};

ChannelUnpremultiplyMapEntry _unpremultiply_i8 ( &channel_unpremultiply_int<int8> );
ChannelUnpremultiplyMapEntry _unpremultiply_u8 ( &channel_unpremultiply_int<uint8> );
ChannelUnpremultiplyMapEntry _unpremultiply_i16( &channel_unpremultiply_int<int16> );
ChannelUnpremultiplyMapEntry _unpremultiply_u16( &channel_unpremultiply_int<uint16> );
ChannelUnpremultiplyMapEntry _unpremultiply_i32( &channel_unpremultiply_int<int32> );
ChannelUnpremultiplyMapEntry _unpremultiply_u32( &channel_unpremultiply_int<uint32> );
ChannelUnpremultiplyMapEntry _unpremultiply_i64( &channel_unpremultiply_int<int64> );
ChannelUnpremultiplyMapEntry _unpremultiply_u64( &channel_unpremultiply_int<uint64> );
ChannelUnpremultiplyMapEntry _unpremultiply_f32( &channel_unpremultiply_float<float> );
ChannelUnpremultiplyMapEntry _unpremultiply_f64( &channel_unpremultiply_float<double> );


//-----------------------------------------------------------------------------------------
// Main conversion functions




void vw::convert( ImageBuffer const& dst, ImageBuffer const& src, bool rescale ) {
  VW_ASSERT( dst.format.cols==src.format.cols && dst.format.rows==src.format.rows,
             ArgumentErr() << "Destination buffer has wrong size." );

  // We only support a few special conversions, and the general case where
  // the source and destination formats are the same.  Below we assume that
  // we're doing a supported conversion, so we check first.
  if( dst.format.pixel_format != src.format.pixel_format ) {
    // We freely convert between multi-channel and multi-plane images,
    // by aliasing the multi-channel buffer as a multi-plane buffer.
    if( src.format.pixel_format==VW_PIXEL_SCALAR && dst.format.planes==1
        && src.format.planes==num_channels( dst.format.pixel_format ) ) {
      ImageBuffer new_dst = dst;
      new_dst.format.pixel_format = VW_PIXEL_SCALAR;
      new_dst.format.planes = src.format.planes;
      new_dst.pstride = channel_size( dst.format.channel_type );
      return convert( new_dst, src );
    }
    else if( dst.format.pixel_format==VW_PIXEL_SCALAR && src.format.planes==1
             && dst.format.planes==num_channels( src.format.pixel_format ) ) {
      ImageBuffer new_src = src;
      new_src.format.pixel_format = VW_PIXEL_SCALAR;
      new_src.format.planes = dst.format.planes;
      new_src.pstride = channel_size( src.format.channel_type );
      return convert( dst, new_src );
    }
    // We support conversions between user specified generic pixel
    // types and the pixel types with an identical number of channels.
    if ( ( src.format.pixel_format == VW_PIXEL_SCALAR_MASKED     && dst.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.format.pixel_format == VW_PIXEL_SCALAR_MASKED     && src.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.format.pixel_format == VW_PIXEL_GRAY_MASKED       && dst.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.format.pixel_format == VW_PIXEL_GRAY_MASKED       && src.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.format.pixel_format == VW_PIXEL_RGB_MASKED        && dst.format.pixel_format == VW_PIXEL_RGBA  ) ||
         ( dst.format.pixel_format == VW_PIXEL_RGB_MASKED        && src.format.pixel_format == VW_PIXEL_RGBA  ) ||
         ( src.format.pixel_format == VW_PIXEL_GENERIC_1_CHANNEL && dst.format.pixel_format == VW_PIXEL_GRAY  ) ||
         ( dst.format.pixel_format == VW_PIXEL_GENERIC_1_CHANNEL && src.format.pixel_format == VW_PIXEL_GRAY  ) ||
         ( src.format.pixel_format == VW_PIXEL_GENERIC_2_CHANNEL && dst.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.format.pixel_format == VW_PIXEL_GENERIC_2_CHANNEL && src.format.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.format.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && dst.format.pixel_format == VW_PIXEL_RGB   ) ||
         ( dst.format.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && src.format.pixel_format == VW_PIXEL_RGB   ) ||
         ( src.format.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && dst.format.pixel_format == VW_PIXEL_XYZ   ) ||
         ( dst.format.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && src.format.pixel_format == VW_PIXEL_XYZ   ) ||
         ( src.format.pixel_format == VW_PIXEL_GENERIC_4_CHANNEL && dst.format.pixel_format == VW_PIXEL_RGBA  ) ||
         ( dst.format.pixel_format == VW_PIXEL_GENERIC_4_CHANNEL && src.format.pixel_format == VW_PIXEL_RGBA  ) ) {
      // Do nothing, these combinations are ok to convert.
    }
    // Other than that, we only support conversion between the core pixel formats
    else if( ( src.format.pixel_format!=VW_PIXEL_GRAY && src.format.pixel_format!=VW_PIXEL_GRAYA &&
               src.format.pixel_format!=VW_PIXEL_RGB  && src.format.pixel_format!=VW_PIXEL_RGBA  &&
               src.format.pixel_format!=VW_PIXEL_XYZ) ||
             ( dst.format.pixel_format!=VW_PIXEL_GRAY && dst.format.pixel_format!=VW_PIXEL_GRAYA &&
               dst.format.pixel_format!=VW_PIXEL_RGB  && dst.format.pixel_format!=VW_PIXEL_RGBA  &&
               dst.format.pixel_format!=VW_PIXEL_XYZ) ) {
      vw_throw( ArgumentErr() << "Source and destination buffers have incompatible pixel formats ("
                << pixel_format_name(src.format.pixel_format) << " vs. " << pixel_format_name(dst.format.pixel_format) << ")." );
    }
  }

  // Gather some stats
  size_t src_channels = num_channels( src.format.pixel_format );
  size_t dst_channels = num_channels( dst.format.pixel_format );
  size_t src_chstride = channel_size( src.format.channel_type );
  size_t dst_chstride = channel_size( dst.format.channel_type );

  int32 copy_length = (src_channels==dst_channels) ? src_channels : (src_channels<3) ? 1 : (dst_channels>=3) ? 3 : 0;

  // Decide how alpha handling and other issues will be done
  bool unpremultiply_src = false, premultiply_src = false, premultiply_dst = false;
  {
    const ImageFormat& srcf = src.format, dstf = dst.format;
    bool src_alpha = (srcf.pixel_format == VW_PIXEL_GRAYA || dstf.pixel_format == VW_PIXEL_RGBA);
    bool dst_alpha = (dstf.pixel_format == VW_PIXEL_GRAYA || dstf.pixel_format == VW_PIXEL_RGBA);
    unpremultiply_src = (src_alpha && srcf.premultiplied && !dstf.premultiplied);
    premultiply_src   = (src_alpha && !dst_alpha && !srcf.premultiplied);
    premultiply_dst   = (src_alpha && dst_alpha && !srcf.premultiplied && dstf.premultiplied);
  }

  bool triplicate = src_channels<3    && dst_channels>=3;
  bool average    = src_channels >=3  && dst_channels<3;
  bool add_alpha  = src_channels%2==1 && dst_channels%2==0;
  bool copy_alpha = src_channels!=dst_channels && src_channels%2==0 && dst_channels%2==0;

  // Get handler functions for all of the input data types.
  // - This could be replaced with function calls containing a switch statement.
  channel_convert_func conv_func = rescale
    ? channel_convert_rescale_map->operator[](std::make_pair(src.format.channel_type,dst.format.channel_type))
    : channel_convert_map->operator[](std::make_pair(src.format.channel_type,dst.format.channel_type));
  channel_set_max_func  max_func = channel_set_max_map->operator[](dst.format.channel_type);
  channel_average_func       avg_func               = channel_average_map->operator[](dst.format.channel_type);
  channel_unpremultiply_func unpremultiply_src_func = channel_unpremultiply_map->operator[](src.format.channel_type);
  channel_premultiply_func   premultiply_src_func   = channel_premultiply_map->operator[](src.format.channel_type);
  channel_premultiply_func   premultiply_dst_func   = channel_premultiply_map->operator[](dst.format.channel_type);
  if( !conv_func || !max_func || !avg_func || !unpremultiply_src_func || !premultiply_dst_func || !premultiply_src_func )
    vw_throw( NoImplErr() << "Unsupported channel type combination in convert (" << src.format.channel_type << ", " << dst.format.channel_type << ")!" );

  int32 max_channels = std::max( src_channels, dst_channels );

  boost::scoped_array<uint8> src_buf(new uint8[max_channels*src_chstride]);
  boost::scoped_array<uint8> dst_buf(new uint8[max_channels*dst_chstride]);

  // Loop through all of the pixels in the source data
  // - Data pointers are always in bytes, will be advanced according to data element size.
  uint8 *src_ptr_p = (uint8*)src.data;
  uint8 *dst_ptr_p = (uint8*)dst.data;
  for( uint32 p=0; p<src.format.planes; ++p ) {
    uint8 *src_ptr_r = src_ptr_p;
    uint8 *dst_ptr_r = dst_ptr_p;
    for( uint32 r=0; r<src.format.rows; ++r ) {
      uint8 *src_ptr_c = src_ptr_r;
      uint8 *dst_ptr_c = dst_ptr_r;
      for( uint32 c=0; c<src.format.cols; ++c ) {

        // Setup the buffers, adjusting premultiplication if needed
        uint8 *src_ptr = src_ptr_c;
        uint8 *dst_ptr = dst_ptr_c;
        if( unpremultiply_src ) {
          unpremultiply_src_func( src_ptr, src_buf.get(), src_channels );
          src_ptr = src_buf.get();
        }
        else if( premultiply_src ) {
          premultiply_src_func( src_ptr, src_buf.get(), src_channels );
          src_ptr = src_buf.get();
        }

        // Copy/convert, unrolling the common multi-channel cases
        if( copy_length==4 ) {
          conv_func( src_ptr,                dst_ptr );
          conv_func( src_ptr+  src_chstride, dst_ptr+  dst_chstride );
          conv_func( src_ptr+2*src_chstride, dst_ptr+2*dst_chstride );
          conv_func( src_ptr+3*src_chstride, dst_ptr+3*dst_chstride );
        }
        else if( copy_length==3 ) {
          conv_func( src_ptr,                dst_ptr );
          conv_func( src_ptr+  src_chstride, dst_ptr+  dst_chstride );
          conv_func( src_ptr+2*src_chstride, dst_ptr+2*dst_chstride );
        }
        else if( copy_length==2) {
          conv_func( src_ptr,              dst_ptr );
          conv_func( src_ptr+src_chstride, dst_ptr+dst_chstride );
        }
        else if( copy_length==1 ) {
          conv_func( src_ptr, dst_ptr );
        }
        else {
          for( int32 ch=0; ch<copy_length; ++ch ) {
            conv_func( src_ptr+ch*src_chstride, dst_ptr+ch*dst_chstride );
          }
        }

        // Handle the special pixel format conversions
        if( triplicate ) {
          // Duplicate the input channel twice more
          conv_func( src_ptr, dst_ptr+  dst_chstride );
          conv_func( src_ptr, dst_ptr+2*dst_chstride );
        }
        else if( average ) {
          for( int32 ch=0; ch<3; ++ch ) {
            conv_func( src_ptr+ch*src_chstride, dst_buf.get()+ch*dst_chstride );
          }
         avg_func( dst_buf.get(), dst_ptr, 3 );
        }
        if( copy_alpha ) {
          conv_func( src_ptr+(src_channels-1)*src_chstride, dst_ptr+(dst_channels-1)*dst_chstride );
        }
        else if( add_alpha ) {
          max_func( dst_ptr+(dst_channels-1)*dst_chstride );
        }

        // Finally, adjust destination premultiplication if needed
        if( premultiply_dst ) {
          premultiply_dst_func( dst_ptr, dst_ptr, dst_channels );
        }
        // Advance pointers by the appropriate size in bytes
        src_ptr_c += src.cstride;
        dst_ptr_c += dst.cstride;
      }
      src_ptr_r += src.rstride;
      dst_ptr_r += dst.rstride;
    }
    src_ptr_p += src.pstride;
    dst_ptr_p += dst.pstride;
  }
} // End function convert


// TODO: Lots of duplicated code here, would be nice to clean up these functions.
void vw::can_convert( ImageFormat const& dst, ImageFormat const& src ) {
  VW_ASSERT( dst.cols==src.cols && dst.rows==src.rows,
             ArgumentErr() << "Destination buffer has wrong size." );

  // We only support a few special conversions, and the general case where
  // the source and destination formats are the same.  Below we assume that
  // we're doing a supported conversion, so we check first.
  if( dst.pixel_format != src.pixel_format ) {
    // We freely convert between multi-channel and multi-plane images,
    // by aliasing the multi-channel buffer as a multi-plane buffer.
    if( src.pixel_format==VW_PIXEL_SCALAR && dst.planes==1
        && src.planes==num_channels( dst.pixel_format ) ) {
      ImageFormat new_dst = dst;
      new_dst.pixel_format = VW_PIXEL_SCALAR;
      new_dst.planes = src.planes;
      return can_convert( new_dst, src );
    }
    else if( dst.pixel_format==VW_PIXEL_SCALAR && src.planes==1
             && dst.planes==num_channels( src.pixel_format ) ) {
      ImageFormat new_src = src;
      new_src.pixel_format = VW_PIXEL_SCALAR;
      new_src.planes = dst.planes;
      return can_convert( dst, new_src );
    }
    // We support conversions between user specified generic pixel
    // types and the pixel types with an identical number of channels.
    if ( ( src.pixel_format == VW_PIXEL_SCALAR_MASKED     && dst.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.pixel_format == VW_PIXEL_SCALAR_MASKED     && src.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.pixel_format == VW_PIXEL_GRAY_MASKED       && dst.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.pixel_format == VW_PIXEL_GRAY_MASKED       && src.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.pixel_format == VW_PIXEL_RGB_MASKED        && dst.pixel_format == VW_PIXEL_RGBA  ) ||
         ( dst.pixel_format == VW_PIXEL_RGB_MASKED        && src.pixel_format == VW_PIXEL_RGBA  ) ||
         ( src.pixel_format == VW_PIXEL_GENERIC_1_CHANNEL && dst.pixel_format == VW_PIXEL_GRAY  ) ||
         ( dst.pixel_format == VW_PIXEL_GENERIC_1_CHANNEL && src.pixel_format == VW_PIXEL_GRAY  ) ||
         ( src.pixel_format == VW_PIXEL_GENERIC_2_CHANNEL && dst.pixel_format == VW_PIXEL_GRAYA ) ||
         ( dst.pixel_format == VW_PIXEL_GENERIC_2_CHANNEL && src.pixel_format == VW_PIXEL_GRAYA ) ||
         ( src.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && dst.pixel_format == VW_PIXEL_RGB   ) ||
         ( dst.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && src.pixel_format == VW_PIXEL_RGB   ) ||
         ( src.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && dst.pixel_format == VW_PIXEL_XYZ   ) ||
         ( dst.pixel_format == VW_PIXEL_GENERIC_3_CHANNEL && src.pixel_format == VW_PIXEL_XYZ   ) ||
         ( src.pixel_format == VW_PIXEL_GENERIC_4_CHANNEL && dst.pixel_format == VW_PIXEL_RGBA  ) ||
         ( dst.pixel_format == VW_PIXEL_GENERIC_4_CHANNEL && src.pixel_format == VW_PIXEL_RGBA  ) ) {
      // Do nothing, these combinations are ok to convert.
    }
    // Other than that, we only support conversion between the core pixel formats
    else if( ( src.pixel_format!=VW_PIXEL_GRAY && src.pixel_format!=VW_PIXEL_GRAYA &&
               src.pixel_format!=VW_PIXEL_RGB  && src.pixel_format!=VW_PIXEL_RGBA  &&
               src.pixel_format!=VW_PIXEL_XYZ) ||
             ( dst.pixel_format!=VW_PIXEL_GRAY && dst.pixel_format!=VW_PIXEL_GRAYA &&
               dst.pixel_format!=VW_PIXEL_RGB  && dst.pixel_format!=VW_PIXEL_RGBA  &&
               dst.pixel_format!=VW_PIXEL_XYZ) ) {
      vw_throw( ArgumentErr() << "Source and destination buffers have incompatible pixel formats ("
                << pixel_format_name(src.pixel_format) << " vs. " << pixel_format_name(dst.pixel_format) << ")." );
    }
  }

  channel_convert_func conv_func = channel_convert_map->operator[](std::make_pair(src.channel_type,dst.channel_type));
  channel_set_max_func max_func  = channel_set_max_map->operator[](dst.channel_type);
  if( !conv_func || !max_func )
    vw_throw( NoImplErr() << "Unsupported channel type combination in convert (" 
                          << src.channel_type << ", " << dst.channel_type << ")!" );

  // If we made it to the end, we can do this conversion.
}

//-----------------------------------------------------------------------------------
// SrcImageResource functions

ImageFormat SrcImageResource::format() const {
  ImageFormat fmt;
  fmt.cols         = this->cols();
  fmt.rows         = this->rows();
  fmt.planes       = this->planes();
  fmt.pixel_format = this->pixel_format();
  fmt.channel_type = this->channel_type();
  return fmt;
}

boost::shared_array<const uint8> SrcImageResource::native_ptr() const {
  boost::shared_array<const uint8> data(new uint8[native_size()]);
  this->read(ImageBuffer(format(), const_cast<uint8*>(data.get())), BBox2i(0,0,cols(),rows()));
  return data;
}

size_t SrcImageResource::native_size() const {
  return channel_size(channel_type()) * num_channels(pixel_format()) * cols() * rows() * planes();
}
