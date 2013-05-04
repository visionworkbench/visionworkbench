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


%module imagemanip
%import "_image.i"

%{
#include <vw/Image.h>
%}

%inline %{
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _transpose( ImageT const& image ) { return transpose(image); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _rotate_180( ImageT const& image ) { return rotate_180( image ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _rotate_90_cw( ImageT const& image ) { return rotate_90_cw( image ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _rotate_90_ccw( ImageT const& image ) { return rotate_90_ccw( image ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _flip_vertical( ImageT const& image ) { return flip_vertical( image ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _flip_horizontal( ImageT const& image ) { return flip_horizontal( image ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _crop( ImageT const& image, int x, int y, int cols, int rows ) { return crop( image, x, y, cols, rows ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _select_col( ImageT const& image, int col ) { return select_col( image, col ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _select_row( ImageT const& image, int row ) { return select_row( image, row ); }
  template <class ImageT> vw::ImageViewRef<typename ImageT::pixel_type> _select_plane( ImageT const& image, int plane ) { return select_plane( image, plane ); }
%}

%define %instantiate_manipulation(cname,ctype,pname,ptype,...)
  %template(transpose) _transpose<vw::ImageView<ptype > >;
  %template(transpose) _transpose<vw::ImageViewRef<ptype > >;
  %template(rotate_180) _rotate_180<vw::ImageView<ptype > >;
  %template(rotate_180) _rotate_180<vw::ImageViewRef<ptype > >;
  %template(rotate_90_cw) _rotate_90_cw<vw::ImageView<ptype > >;
  %template(rotate_90_cw) _rotate_90_cw<vw::ImageViewRef<ptype > >;
  %template(rotate_90_ccw) _rotate_90_ccw<vw::ImageView<ptype > >;
  %template(rotate_90_ccw) _rotate_90_ccw<vw::ImageViewRef<ptype > >;
  %template(flip_vertical) _flip_vertical<vw::ImageView<ptype > >;
  %template(flip_vertical) _flip_vertical<vw::ImageViewRef<ptype > >;
  %template(flip_horizontal) _flip_horizontal<vw::ImageView<ptype > >;
  %template(flip_horizontal) _flip_horizontal<vw::ImageViewRef<ptype > >;
  %template(crop) _crop<vw::ImageView<ptype > >;
  %template(crop) _crop<vw::ImageViewRef<ptype > >;
  %template(select_col) _select_col<vw::ImageView<ptype > >;
  %template(select_col) _select_col<vw::ImageViewRef<ptype > >;
  %template(select_row) _select_row<vw::ImageView<ptype > >;
  %template(select_row) _select_row<vw::ImageViewRef<ptype > >;
  %template(select_plane) _select_plane<vw::ImageView<ptype > >;
  %template(select_plane) _select_plane<vw::ImageViewRef<ptype > >;
%enddef

%instantiate_for_pixel_types(instantiate_manipulation)
