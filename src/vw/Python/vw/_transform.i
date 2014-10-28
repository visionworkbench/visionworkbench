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


%module(directors="1") transform
%include "std_string.i"
%include "std_vector.i"
%import "_pixel.i"
%import "_image.i"

%{
#define SWIG_FILE_WITH_INIT
#include <vw/Image.h>
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}

namespace vw {
  class NearestPixelInterpolation {};
  class BilinearInterpolation {};
  class BicubicInterpolation {};
}

%define %instantiate_for_pixel_types_and_interpolations(macro,args...)
  %instantiate_for_pixel_types(macro,vw::NearestPixelInterpolation,args)
  %instantiate_for_pixel_types(macro,vw::BilinearInterpolation,args)
  %instantiate_for_pixel_types(macro,vw::BicubicInterpolation,args)
%enddef

%typemap(typecheck) (double* IN_ARRAY2, int DIM1, int DIM2) {
  $1 = is_array($input);
}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};

%exception _homography {
  try {
    $action
  } catch (const vw::ArgumentErr& e) {
    PyErr_SetString(PyExc_TypeError, const_cast<char*>(e.what()));
    SWIG_fail;
  }
}

// These functions assume that the source ImageViewRef has already been edge extended!
%inline %{
  template <class PixelT, class InterpT>
  vw::TransformView<vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT>, vw::TransformRef>
  inline __transform( vw::ImageViewRef<PixelT> const& image, vw::TransformRef const& transform, vw::int32 x_size, vw::int32 y_size, InterpT const& interp ) {
    vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT> interpolated( image, interp );
    return vw::TransformView<vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT>, vw::TransformRef>( interpolated, transform, x_size, y_size );
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _resample( vw::ImageViewRef<PixelT> const& image, double x_scale, double y_scale, InterpT const& interp ) {
    vw::TransformRef txform( vw::ResampleTransform(x_scale, y_scale) );
    return __transform( image, txform, int(image.cols()*x_scale+0.5), int(image.rows()*y_scale+0.5), interp );
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _resize( vw::ImageViewRef<PixelT> const& image, vw::int32 x_size, vw::int32 y_size, InterpT const& interp ) {
    vw::TransformRef txform( vw::ResampleTransform(x_size/(double)image.cols(), y_size/(double)image.rows()) );
    return __transform( image, txform, x_size, y_size, interp );
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _translate( vw::ImageViewRef<PixelT> const& image, double x, double y, InterpT const& interp ) {
    vw::TransformRef txform( vw::TranslateTransform(x,y) );
    return __transform( image, txform, image.cols(), image.rows(), interp );
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _rotate( vw::ImageViewRef<PixelT> const& image, double theta, InterpT const& interp ) {
    vw::RotateTransform rot(theta, vw::Vector2(0,0));
    vw::Vector2 center( (image.cols()-1)/2.0, (image.rows()-1)/2.0 );
    vw::Vector2 offset = center - rot.forward(center);
    vw::TransformRef txform( compose(vw::TranslateTransform(offset.x(),offset.y()),rot) );
    return __transform( image, txform, image.cols(), image.rows(), interp );
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _homography( vw::ImageViewRef<PixelT> const& image, double* matrix, int rows, int cols, vw::int32 x_size, vw::int32 y_size, InterpT const& interp ) {
    if( rows != 3 || cols != 3 ) {
      throw vw::ArgumentErr() << "Array must have dimensions 3x3.  Given array has dimensions " << rows << "x" << cols << ".";
    }
    vw::Matrix<double,3,3> homog_matrix( matrix );
    vw::HomographyTransform homog( homog_matrix );
    vw::TransformRef txform( homog );
    return __transform( image, txform, x_size, y_size, interp );
  }

  //template <class PixelT, class InterpT>
  //vw::ImageViewRef<PixelT> _transform( vw::ImageViewRef<PixelT> const& image, VirtualTransformBase* txform, vw::int32 x_size, vw::int32 y_size, InterpT const& interp ) {
  //
  //}
%}

%define %instantiate_transforms(cname,ctype,pname,ptype,in,...)
  %template(_resample) _resample<ptype, in>;
  %template(_resize) _resize<ptype, in>;
  %template(_translate) _translate<ptype, in>;
  %template(_rotate) _rotate<ptype, in>;
  %template(_homography) _homography<ptype, in>;
%enddef

%instantiate_for_pixel_types_and_interpolations(instantiate_transforms)

%pythoncode {
  from image import edge_extend, ZeroEdgeExtension, ConstantEdgeExtension

  def resample(image,xfactor,yfactor=None,edge=None,interp=None):
    if yfactor is None: yfactor = xfactor
    if edge is None: edge = ConstantEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _resample(edge_extend(image,edge=edge),xfactor,yfactor,interp)

  def resize(image,xsize,ysize,edge=None,interp=None):
    if edge is None: edge = ConstantEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _resize(edge_extend(image,edge=edge),xsize,ysize,interp)

  def translate(image,xoff,yoff,edge=None,interp=None):
    if edge is None: edge = ZeroEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _translate(edge_extend(image,edge=edge),xoff,yoff,interp)

  def rotate(image,theta,edge=None,interp=None):
    if edge is None: edge = ZeroEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _rotate(edge_extend(image,edge=edge),theta,interp)

  def homography(image,matrix,xsize=None,ysize=None,edge=None,interp=None):
    if xsize is None: xsize = image.cols
    if ysize is None: ysize = image.rows
    if edge is None: edge = ZeroEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _homography(edge_extend(image,edge=edge),matrix,xsize,ysize,interp)
}
