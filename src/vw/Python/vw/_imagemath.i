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


%module imagemath
%import "_image.i"

%{
#include <vw/Image.h>
%}

%inline %{
  template <class PixelT> vw::ImageViewRef<PixelT> _add_image_image( vw::ImageViewRef<PixelT> const& image1, vw::ImageViewRef<PixelT> const& image2 ) { return image1 + image2; }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _add_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return image + value; }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _radd_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return value + image; }
  template <class PixelT> vw::ImageView<PixelT> _iadd_image_image( vw::ImageView<PixelT> const& image1, vw::ImageViewRef<PixelT> const& image2 ) { return image1 += image2; }
  template <class PixelT, class ValueT> vw::ImageView<PixelT> _iadd_image_value( vw::ImageView<PixelT> const& image, ValueT value ) { return image += value; }

  template <class PixelT> vw::ImageViewRef<PixelT> _sub_image_image( vw::ImageViewRef<PixelT> const& image1, vw::ImageViewRef<PixelT> const& image2 ) { return image1 - image2; }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _sub_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return image - value; }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _rsub_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return value - image; }
  template <class PixelT> vw::ImageView<PixelT> _isub_image_image( vw::ImageView<PixelT> const& image1, vw::ImageViewRef<PixelT> const& image2 ) { return image1 -= image2; }
  template <class PixelT, class ValueT> vw::ImageView<PixelT> _isub_image_value( vw::ImageView<PixelT> const& image, ValueT value ) { return image -= value; }

  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _mul_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return vw::pixel_cast<PixelT>( image * value ); }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _rmul_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return vw::pixel_cast<PixelT>( value * image ); }
  template <class PixelT, class ValueT> vw::ImageView<PixelT> _imul_image_value( vw::ImageView<PixelT> const& image, ValueT value ) { return image *= value; }

  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _div_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return vw::pixel_cast<PixelT>( image / value ); }
  template <class PixelT, class ValueT> vw::ImageViewRef<PixelT> _rdiv_image_value( vw::ImageViewRef<PixelT> const& image, ValueT value ) { return vw::pixel_cast<PixelT>( value / image ); }
  template <class PixelT, class ValueT> vw::ImageView<PixelT> _idiv_image_value( vw::ImageView<PixelT> const& image, ValueT value ) { return image /= value; }

  // These functions are only designed to be instantiated on floating-point pixel types.
  template <class PixelT> vw::ImageViewRef<PixelT> _sin( vw::ImageViewRef<PixelT> const& image ) { return vw::sin(image); }
  template <class PixelT> vw::ImageViewRef<PixelT> _cos( vw::ImageViewRef<PixelT> const& image ) { return vw::cos(image); }
  template <class PixelT> vw::ImageViewRef<PixelT> _tan( vw::ImageViewRef<PixelT> const& image ) { return vw::tan(image); }
%}

%define %instantiate_math_unary_func(pname,ptype,fname)
  %template(_##fname) _##fname<ptype >;
  %pythoncode {
    def fname(image):
      return _##fname(image.ref())
  }
%enddef

%define %instantiate_math(cname,ctype,pname,ptype,...)
#if %is_float(ctype)
  %instantiate_math_unary_func(pname,ptype,sin)
  %instantiate_math_unary_func(pname,ptype,cos)
  %instantiate_math_unary_func(pname,ptype,tan)
#endif
  %template(image_##pname##_add) _add_image_image<ptype >;
  %template(image_##pname##_add) _add_image_value<ptype, ptype >;
  %template(image_##pname##_add) _add_image_value<ptype, ctype >;
  %template(image_##pname##_radd) _radd_image_value<ptype, ptype >;
  %template(image_##pname##_radd) _radd_image_value<ptype, ctype >;
  %template(image_##pname##_iadd) _iadd_image_image<ptype >;
  %template(image_##pname##_iadd) _iadd_image_value<ptype, ptype >;
  %template(image_##pname##_iadd) _iadd_image_value<ptype, ctype >;

  %template(image_##pname##_sub) _sub_image_image<ptype >;
  %template(image_##pname##_sub) _sub_image_value<ptype, ptype >;
  %template(image_##pname##_sub) _sub_image_value<ptype, ctype >;
  %template(image_##pname##_rsub) _rsub_image_value<ptype, ptype >;
  %template(image_##pname##_rsub) _rsub_image_value<ptype, ctype >;
  %template(image_##pname##_isub) _isub_image_image<ptype >;
  %template(image_##pname##_isub) _isub_image_value<ptype, ptype >;
  %template(image_##pname##_isub) _isub_image_value<ptype, ctype >;

  %template(image_##pname##_mul) _mul_image_value<ptype, vw::float32 >;
  %template(image_##pname##_rmul) _rmul_image_value<ptype, vw::float32 >;
  %template(image_##pname##_imul) _imul_image_value<ptype, vw::float32 >;

  %template(image_##pname##_div) _div_image_value<ptype, vw::float32 >;
  %template(image_##pname##_rdiv) _rdiv_image_value<ptype, vw::float32 >;
  %template(image_##pname##_idiv) _idiv_image_value<ptype, vw::float32 >;

  %pythoncode {
    image.ImageView_##pname.__add__ = lambda self, other: image_##pname##_add( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__add__ = lambda self, other: image_##pname##_add( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__radd__ = lambda self, other: image_##pname##_radd( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__radd__ = lambda self, other: image_##pname##_radd( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__iadd__ = lambda self, other: image_##pname##_iadd( self, image._ref_if_image(other) )

    image.ImageView_##pname.__sub__ = lambda self, other: image_##pname##_sub( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__sub__ = lambda self, other: image_##pname##_sub( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__rsub__ = lambda self, other: image_##pname##_rsub( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__rsub__ = lambda self, other: image_##pname##_rsub( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__isub__ = lambda self, other: image_##pname##_isub( self, image._ref_if_image(other) )

    image.ImageView_##pname.__mul__ = lambda self, other: image_##pname##_mul( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__mul__ = lambda self, other: image_##pname##_mul( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__rmul__ = lambda self, other: image_##pname##_rmul( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__rmul__ = lambda self, other: image_##pname##_rmul( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__imul__ = lambda self, other: image_##pname##_imul( self, image._ref_if_image(other) )

    image.ImageView_##pname.__rdiv__ = lambda self, other: image_##pname##_rdiv( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__rdiv__ = lambda self, other: image_##pname##_rdiv( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__rdiv__ = lambda self, other: image_##pname##_rdiv( self.ref(), image._ref_if_image(other) )
    image.ImageViewRef_##pname.__rdiv__ = lambda self, other: image_##pname##_rdiv( self.ref(), image._ref_if_image(other) )
    image.ImageView_##pname.__idiv__ = lambda self, other: image_##pname##_idiv( self, image._ref_if_image(other) )
  }
%enddef

%instantiate_for_pixel_types(instantiate_math)
