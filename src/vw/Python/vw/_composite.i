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


%module composite
%import "_image.i"

%{
#include <vw/Image.h>
#include <vw/Mosaic.h>
%}

namespace vw {
namespace mosaic {

  template <class PixelT>
  class ImageComposite {
  public:
    ImageComposite();


    void prepare();
    void prepare( vw::BBox2i const& total_bbox );

    void set_draft_mode( bool draft_mode );

    %extend {
      int get_cols() const { return self->cols(); }
      int get_rows() const { return self->rows(); }
      int get_planes() const { return self->planes(); }
      int get_channels() const { return self->channels(); }
      vw::BBox2i const& get_bbox() const { return self->bbox(); }
      vw::BBox2i const& get_source_data_bbox() const { return self->source_data_bbox(); }
      vw::ImageViewRef<PixelT> ref() const { return *self; }

      void _insert( vw::ImageViewRef<PixelT> const& image, int x, int y ) { self->insert(image,x,y); }
    }

    %pythoncode {
      cols = property(get_cols)
      rows = property(get_rows)
      planes = property(get_planes)
      channels = property(get_channels)
      bbox = property(get_bbox)
      source_data_bbox = property(get_source_data_bbox)
      draft_mode = property(fset=set_draft_mode)

      def insert(self, image, x, y):
        self._insert(image.ref(),x,y)
    }

  };

} // namespace mosaic
} // namespace vw

%pythoncode {

  class ImageComposite(object):

    _pixel_type_table = dict()

    def __init__(self, ptype=None, pformat=None, ctype=None):
      ptype = pixel._compute_pixel_type(pixel.PixelRGBA_float32,ptype,pformat,ctype)
      self.__dict__['_delegate'] = ImageComposite._pixel_type_table[ptype]()

    def __getattr__(self,name):
      return getattr(self._delegate,name)

    def __setattr__(self,name,value):
      return setattr(self._delegate,name,value)
}

%define %instantiate_imagecomposite_types(cname,ctype,pname,ptype,...)
  %template(ImageComposite_##pname) vw::mosaic::ImageComposite<ptype >;
  %pythoncode {
    ImageComposite._pixel_type_table[pixel.pname] = ImageComposite_##pname
    ImageComposite_##pname.pixel_type = pixel.pname
    ImageComposite_##pname.channel_type = pixel.cname
  }
%enddef

%instantiate_for_pixel_types(instantiate_imagecomposite_types)
