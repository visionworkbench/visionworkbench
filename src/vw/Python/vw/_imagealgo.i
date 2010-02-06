// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


%module imagealgo
%import "_image.i"

%{
#include <vw/Image.h>
%}

%inline %{
  template <class ImageT> void _fill( ImageT const& image, typename ImageT::pixel_type value ) {
    return fill( image, value );
  }
  template <class ImageT, class ChannelT> vw::ImageViewRef<typename ImageT::pixel_type> _clamp( ImageT const& image, ChannelT low, ChannelT high ) {
    return clamp( image, low, high );
  }
  template <class ImageT, class ChannelT> vw::ImageViewRef<typename ImageT::pixel_type> _normalize( ImageT const& image, ChannelT low, ChannelT high ) {
    return normalize( image, low, high );
  }
  template <class ImageT, class ChannelT> vw::ImageViewRef<typename ImageT::pixel_type> _threshold( ImageT const& image, ChannelT thresh, ChannelT low, ChannelT high ) {
    return threshold( image, thresh, low, high );
  }
  template <class ImageT> bool _is_opaque( ImageT const& image ) {
    return is_opaque( image );
  }
%}

%define %instantiate_algorithms(cname,ctype,pname,ptype,...)
  %template(fill) _fill<vw::ImageView<ptype > >;
  %template(_clamp) _clamp<vw::ImageView<ptype >, ctype >;
  %template(_clamp) _clamp<vw::ImageViewRef<ptype >, ctype >;
  %template(_normalize) _normalize<vw::ImageView<ptype >, ctype >;
  %template(_normalize) _normalize<vw::ImageViewRef<ptype >, ctype >;
  %template(_threshold) _threshold<vw::ImageView<ptype >, ctype >;
  %template(_threshold) _threshold<vw::ImageViewRef<ptype >, ctype >;
  %template(_is_opaque) _is_opaque<vw::ImageView<ptype > >;
  %template(_is_opaque) _is_opaque<vw::ImageViewRef<ptype > >;
%enddef

%instantiate_for_pixel_types(instantiate_algorithms)

%pythoncode {
  def fill( image ):
    _fill( image )

  def clamp( image, low=None, high=None ):
    cmin, cmax = pixel.channel_range( image.pixel_type )
    return _clamp( image, low or cmin, high or cmax )

  def normalize( image, low=None, high=None ):
    cmin, cmax = pixel.channel_range( image.pixel_type )
    return _normalize( image, low or cmin, high or cmax )

  def threshold( image, thresh=None, low=None, high=None ):
    cmin, cmax = pixel.channel_range( image.pixel_type )
    return _threshold( image, thresh or cmin, low or cmin, high or cmax )

  def is_opaque( image ):
    return _is_opaque( image )
}
