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


%module pixelcast
%import "_image.i"

%{
#include <vw/Image.h>
%}

%inline %{

  template <class OutPixelT, class InImageT>
  vw::ImageViewRef<OutPixelT> _pixel_cast( InImageT const& image ) {
    return vw::pixel_cast<OutPixelT>( image );
  }

%}

%pythoncode {

  _pixel_cast_table = dict()

  def pixel_cast(object, ptype=None, pformat=None, ctype=None):
    ptype = pixel._compute_pixel_type(None,ptype,pformat,ctype)
    return _pixel_cast_table[ptype](object)

}

%define %instantiate_pixelcast_helper(cname1,ctype1,cname2,ctype2,pname2,ptype2,args...)
  %template(_pixel_cast_##pname2##_##cname2) _pixel_cast<ptype2<ctype2>, vw::ImageView<ptype2<ctype1> > >;
  %template(_pixel_cast_##pname2##_##cname2) _pixel_cast<ptype2<ctype2>, vw::ImageViewRef<ptype2<ctype1> > >;
%enddef

%define %instantiate_pixelcast(cname,ctype,pname,ptype,args...)
  %instantiate_for_channel_types(instantiate_pixelcast_helper,cname,ctype,pname,ptype)
  %pythoncode {
    _pixel_cast_table[pixel.pname##_##cname] = _pixelcast._pixel_cast_##pname##_##cname
  }
%enddef

#define %instantiate_pixelcast_for_pixel_format(pname,ptype,args...) %instantiate_for_channel_types(instantiate_pixelcast,pname,ptype)
%instantiate_for_pixel_formats( instantiate_pixelcast_for_pixel_format )

