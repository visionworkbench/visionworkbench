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


%module filter
%include "std_string.i"
%include "std_vector.i"
%import "_pixel.i"
%import "_image.i"

%{
#include <vw/Image.h>
%}

%inline %{
  template <class PixelT, class KernelT, class EdgeT>
  vw::ImageViewRef<PixelT> _convolution_filter( vw::ImageViewRef<PixelT> const& image, vw::ImageView<KernelT> const& kernel, vw::int32 cx, vw::int32 cy, EdgeT const& edge ) {
    return vw::pixel_cast<PixelT>( convolution_filter( image, kernel, cx, cy, edge ) );
  }
  template <class PixelT, class KernelT, class EdgeT>
  vw::ImageViewRef<PixelT> _separable_convolution_filter( vw::ImageViewRef<PixelT> const& image, std::vector<KernelT> const& xkernel, std::vector<KernelT> const& ykernel, vw::int32 cx, vw::int32 cy, EdgeT const& edge ) {
    return vw::pixel_cast<PixelT>( separable_convolution_filter( image, xkernel, ykernel, cx, cy, edge ) );
  }
  template <class PixelT, class EdgeT>
  vw::ImageViewRef<PixelT> _gaussian_filter( vw::ImageViewRef<PixelT> const& image, double xsigma, double ysigma, vw::int32 xkernel, vw::int32 ykernel, EdgeT const& edge ) {
    return vw::pixel_cast<PixelT>( gaussian_filter( image, xsigma, ysigma, xkernel, ykernel, edge ) );
  }
  template <class PixelT, class EdgeT>
  vw::ImageViewRef<PixelT> _derivative_filter( vw::ImageViewRef<PixelT> const& image, vw::int32 xderiv, vw::int32 yderiv, vw::int32 xkernel, vw::int32 ykernel, EdgeT const& edge ) {
    return vw::pixel_cast<PixelT>( derivative_filter( image, xderiv, yderiv, xkernel, ykernel, edge ) );
  }
  template <class PixelT, class EdgeT>
  vw::ImageViewRef<PixelT> _laplacian_filter( vw::ImageViewRef<PixelT> const& image, EdgeT const& edge ) {
    return vw::pixel_cast<PixelT>( laplacian_filter( image, edge ) );
  }
%}

%define %instantiate_filters(cname,ctype,pname,ptype,ee,...)
  %template(_convolution_filter) _convolution_filter<ptype, vw::float32, ee>;
  %template(_separable_convolution_filter) _separable_convolution_filter<ptype, vw::float32, ee>;
  %template(_gaussian_filter) _gaussian_filter<ptype, ee>;
  %template(_derivative_filter) _derivative_filter<ptype, ee>;
  %template(_laplacian_filter) _laplacian_filter<ptype, ee>;
%enddef

%instantiate_for_pixel_types_and_edge_extensions(instantiate_filters)

%pythoncode {
  from image import ConstantEdgeExtension

  def convolution_filter(image,kernel,cx=None,cy=None,edge=None):
    if cx is None: cx = int( kernel.cols / 2 )
    if cy is None: cy = int( kernel.rows / 2 )
    if edge is None: edge = ConstantEdgeExtension()
    return _convolution_filter(image.ref(),kernel,cx,cy,edge)

  def separable_convolution_filter(image,xkernel=(),ykernel=(),cx=None,cy=None,edge=None):
    if len(xkernel) == 0 and len(ykernel) == 0: return image
    if cx is None: cx = int( len(xkernel) / 2 )
    if cy is None: cy = int( len(ykernel) / 2 )
    if edge is None: edge = ConstantEdgeExtension()
    return _separable_convolution_filter(image.ref(),xkernel,ykernel,cx,cy,edge)

  def gaussian_filter(image,xsigma,ysigma=None,xsize=None,ysize=None,edge=None):
    if ysigma is None: ysigma = xsigma
    if xsize is None: xsize = 0
    if ysize is None: ysize = xsize
    if edge is None: edge = ConstantEdgeExtension()
    return _gaussian_filter(image.ref(),xsigma,ysigma,xsize,ysize,edge)

  def derivative_filter(image,xderiv,yderiv,xsize=None,ysize=None,edge=None):
    if xsize is None: xsize = 0
    if ysize is None: ysize = 0
    if edge is None: edge = ConstantEdgeExtension()
    return _derivative_filter(image.ref(),xderiv,yderiv,xsize,ysize,edge)

  def laplacian_filter(image,edge=None):
    if edge is None: edge = ConstantEdgeExtension()
    return _laplacian_filter(image,edge)
}
