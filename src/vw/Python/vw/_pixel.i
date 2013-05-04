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


%module pixel
%include "std_string.i"

%{
#include <vw/Image/PixelTypes.h>
%}

%define %instantiate_for_channel_types(macro,args...)
  %macro( uint8,   vw::uint8,   args )
  %macro( int16,   vw::int16,   args )
  %macro( uint16,  vw::uint16,  args )
  %macro( float32, vw::float32, args )
%enddef

%define %instantiate_for_pixel_formats(macro,args...)
  %macro( PixelGray,  vw::PixelGray,  args )
  %macro( PixelGrayA, vw::PixelGrayA, args )
  %macro( PixelRGB,   vw::PixelRGB,   args )
  %macro( PixelRGBA,  vw::PixelRGBA,  args )
%enddef

#define %instantiate_for_scalar_pixel_for_channel_type(cname,ctype,macro,args...) %macro(cname,ctype,cname,ctype,args)
#define %instantiate_for_pixel_type_for_channel_type(cname,ctype,macro,pname,ptype,args...) %macro(cname,ctype,pname##_##cname,ptype<ctype>,args)
#define %instantiate_for_pixel_type_for_channel_types(pname,ptype,macro,args...) %instantiate_for_channel_types(instantiate_for_pixel_type_for_channel_type,macro,pname,ptype,args)

%define %instantiate_for_pixel_types(macro,args...)
  %instantiate_for_channel_types(instantiate_for_scalar_pixel_for_channel_type,macro,args)
  %instantiate_for_pixel_formats(instantiate_for_pixel_type_for_channel_types,macro,args)
%enddef

#define %is_float(ctype) #ctype=="vw::float16" || #ctype=="vw::float32" || #ctype=="vw::float64"

// This duplication is annoying.
namespace vw {
  enum PixelFormatEnum {
    VW_PIXEL_UNKNOWN = 0,
    VW_PIXEL_SCALAR = 1,
    VW_PIXEL_GRAY = 2,
    VW_PIXEL_GRAYA = 3,
    VW_PIXEL_RGB = 4,
    VW_PIXEL_RGBA = 5,
    VW_PIXEL_HSV = 6,
    VW_PIXEL_XYZ = 7,
    VW_PIXEL_USER = 100
  };

  enum ChannelTypeEnum {
    VW_CHANNEL_UNKNOWN = 0,
    VW_CHANNEL_INT8 = 1,
    VW_CHANNEL_UINT8 = 2,
    VW_CHANNEL_INT16 = 3,
    VW_CHANNEL_UINT16 = 4,
    VW_CHANNEL_INT32 = 5,
    VW_CHANNEL_UINT32 = 6,
    VW_CHANNEL_INT64 = 7,
    VW_CHANNEL_UINT64 = 8,
    VW_CHANNEL_FLOAT16 = 9,
    VW_CHANNEL_FLOAT32 = 10,
    VW_CHANNEL_FLOAT64 = 11,
    VW_CHANNEL_BOOL = 12,
    VW_CHANNEL_CHAR = 13,
    VW_CHANNEL_USER = 100
  };

  typedef signed char int8;
  typedef unsigned char uint8;
  typedef signed short int16;
  typedef unsigned short uint16;
  typedef signed int int32;
  typedef unsigned int uint32;
  typedef signed long long int64;
  typedef unsigned long long uint64;
  typedef float float32;
  typedef double float64;
}

%pythoncode{
  _channel_type_table = {}
  _channel_range_table = {}

  def channel_type(object):
    try:
      return object.channel_type
    except:
      return object

  def channel_range(object):
    return _channel_range_table[ channel_type(object) ]

  def pixel_type(object):
    return object.pixel_type

  def pixel_format(object):
    return pixel_type(object).pixel_format

  class PixelFormat(dict):
    def __init__(self, name):
      dict.__init__(self)
      self.name = name
    def __repr__(self):
      return "<pixel format 'vw.%s'>" % self.name

  PixelScalar = PixelFormat("Scalar")

  _pixel_format_table = { VW_PIXEL_SCALAR : PixelScalar }
}

%define %instantiate_channel_type(cname,cid,crange)
  %pythoncode {
    from numpy import cname
    _channel_type_table[cid] = cname
    _channel_range_table[cname] = crange
    PixelScalar[cname] = cname
  }
%enddef

%instantiate_channel_type(uint8, VW_CHANNEL_UINT8, (0,255))
%instantiate_channel_type(int16, VW_CHANNEL_INT16, (0,32767))
%instantiate_channel_type(uint16, VW_CHANNEL_UINT16, (0,65535))
%instantiate_channel_type(float32, VW_CHANNEL_FLOAT32, (0.0,1.0))

%define %instantiate_pixel_type_for_channel_type(cname,ctype,pname,ptype)
  %template(pname##_##cname) ptype<ctype >;
  %pythoncode {
    pname##_##cname.channel_type = cname
    pname##_##cname.pixel_format = pname
    pname[cname] = pname##_##cname
  }
%enddef

%define %instantiate_pixel_type(pname,ptype,pid)
  %pythoncode {
    pname = PixelFormat(#pname)
    _pixel_format_table[pid] = pname
  }
  %instantiate_for_channel_types(instantiate_pixel_type_for_channel_type,pname,ptype)
%enddef

namespace vw {
  template <class ChannelT>
  class PixelGray {
  public:
    PixelGray();
    PixelGray( ChannelT v );
    %extend {
      ChannelT get_v() { return self->v(); }
      void set_v( ChannelT v ) { self->v() = v; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if name is 'v': return self.get_v()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if name is 'v': self.set_v(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelGray,vw::PixelGray,VW_PIXEL_GRAY)

namespace vw {
  template <class ChannelT>
  class PixelGrayA {
  public:
    PixelGrayA();
    PixelGrayA( ChannelT v );
    PixelGrayA( ChannelT v, ChannelT a );
    %extend {
      ChannelT get_v() { return self->v(); }
      void set_v( ChannelT v ) { self->v() = v; }
      ChannelT get_a() { return self->a(); }
      void set_a( ChannelT a ) { self->a() = a; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if   name is 'v': return self.get_v()
        elif name is 'a': return self.get_a()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if   name is 'v': self.set_v(value)
        elif name is 'a': self.set_a(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelGrayA,vw::PixelGrayA,VW_PIXEL_GRAYA)

namespace vw {
  template <class ChannelT>
  class PixelRGB {
  public:
    PixelRGB();
    PixelRGB( ChannelT v );
    PixelRGB( ChannelT r, ChannelT g, ChannelT b );
    %extend {
      ChannelT get_r() { return self->r(); }
      void set_r( ChannelT val ) { self->r() = val; }
      ChannelT get_g() { return self->g(); }
      void set_g( ChannelT val ) { self->g() = val; }
      ChannelT get_b() { return self->b(); }
      void set_b( ChannelT val ) { self->b() = val; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if   name is 'r': return self.get_r()
        elif name is 'g': return self.get_g()
        elif name is 'b': return self.get_b()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if   name is 'r': self.set_r(value)
        elif name is 'g': self.set_g(value)
        elif name is 'b': self.set_b(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelRGB,vw::PixelRGB,VW_PIXEL_RGB)

namespace vw {
  template <class ChannelT>
  class PixelRGBA {
  public:
    PixelRGBA();
    PixelRGBA( ChannelT v );
    PixelRGBA( ChannelT r, ChannelT g, ChannelT b, ChannelT a );
    %extend {
      ChannelT get_r() { return self->r(); }
      void set_r( ChannelT val ) { self->r() = val; }
      ChannelT get_g() { return self->g(); }
      void set_g( ChannelT val ) { self->g() = val; }
      ChannelT get_b() { return self->b(); }
      void set_b( ChannelT val ) { self->b() = val; }
      ChannelT get_a() { return self->a(); }
      void set_a( ChannelT val ) { self->a() = val; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if   name is 'r': return self.get_r()
        elif name is 'g': return self.get_g()
        elif name is 'b': return self.get_b()
        elif name is 'a': return self.get_a()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if   name is 'r': self.set_r(value)
        elif name is 'g': self.set_g(value)
        elif name is 'b': self.set_b(value)
        elif name is 'a': self.set_a(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelRGBA,vw::PixelRGBA,VW_PIXEL_RGBA)

namespace vw {
  template <class ChannelT>
  class PixelHSV {
  public:
    PixelHSV();
    PixelHSV( ChannelT v );
    PixelHSV( ChannelT h, ChannelT s, ChannelT v );
    %extend {
      ChannelT get_h() { return self->h(); }
      void set_h( ChannelT val ) { self->h() = val; }
      ChannelT get_s() { return self->s(); }
      void set_s( ChannelT val ) { self->s() = val; }
      ChannelT get_v() { return self->v(); }
      void set_v( ChannelT val ) { self->v() = val; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if   name is 'h': return self.get_h()
        elif name is 's': return self.get_s()
        elif name is 'v': return self.get_v()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if   name is 'h': self.set_h(value)
        elif name is 's': self.set_s(value)
        elif name is 'v': self.set_v(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelHSV,vw::PixelHSV,VW_PIXEL_HSV)

namespace vw {
  template <class ChannelT>
  class PixelXYZ {
  public:
    PixelXYZ();
    PixelXYZ( ChannelT v );
    PixelXYZ( ChannelT x, ChannelT y, ChannelT z );
    %extend {
      ChannelT get_x() { return self->x(); }
      void set_x( ChannelT val ) { self->x() = val; }
      ChannelT get_y() { return self->y(); }
      void set_y( ChannelT val ) { self->y() = val; }
      ChannelT get_z() { return self->z(); }
      void set_z( ChannelT val ) { self->z() = val; }
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
    %pythoncode{
      _old_getattr = __getattr__
      def __getattr__(self,name):
        if   name is 'x': return self.get_x()
        elif name is 'y': return self.get_y()
        elif name is 'z': return self.get_z()
        else: return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if   name is 'x': self.set_x(value)
        elif name is 'y': self.set_y(value)
        elif name is 'z': self.set_z(value)
        else: self._old_setattr(name,value)
    %}
  };
}
%instantiate_pixel_type(PixelXYZ,vw::PixelXYZ,VW_PIXEL_XYZ)

%pythoncode {
  def _compute_pixel_type(default=None,ptype=None,pformat=None,ctype=None):
    if ptype is None:
      if ctype is None:
        if default is None:
          raise Exception, "No channel type specified"
        ctype = default.channel_type
      if pformat is None:
        if default is None:
          raise Exception, "No pixel format specified"
        pformat = default.pixel_format
      ptype = pformat[ctype]
    elif pformat is not None or ctype is not None:
      raise Exception, "Cannot specify both ptype and pformat/ctype"
    return ptype
}
