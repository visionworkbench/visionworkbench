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


%module image
%include "std_string.i"
%include "std_vector.i"
%import "_pixel.i"

%{
#include <vw/Image.h>
%}

%template(vector_float32) std::vector<vw::float32>;

%pythoncode {
  def isimage(im):
    try:
      im.pixel_type
      return True
    except:
      return False

  def _ref_if_image(im):
    if isimage(im): return im.ref()
    else: return im

  class col(object):
    def __init__(self,i):
      self.index = i;

  class row(object):
    def __init__(self,i):
      self.index = i;

  class plane(object):
    def __init__(self,i):
      self.index = i;

  def Image_getitem(image, pos):
    if isinstance(pos, slice):
      return image.get_region(pos.start[0], pos.start[1], pos.stop[0]-pos.start[0], pos.stop[1]-pos.start[1])
    elif isinstance(pos, col):
      return image.get_col(pos.index)
    elif isinstance(pos, row):
      return image.get_row(pos.index)
    elif isinstance(pos, plane):
      return image.get_plane(pos.index)
    else:
      return image.get_pixel(*pos)

  def Image_setitem(image, pos, val):
    if isinstance(pos, slice):
      image.set_region(val, pos.start[0], pos.start[1], pos.stop[0]-pos.start[0], pos.stop[1]-pos.start[1])
    elif isinstance(pos,col):
      image.set_col(val, pos.index)
    elif isinstance(pos,row):
      image.set_row(val, pos.index)
    elif isinstance(pos,plane):
      image.set_plane(val, pos.index)
    else:
      image.set_pixel(val, *pos)
}

namespace vw {
  class ImageFormat {
  public:
    int32 cols, rows, planes;
    PixelFormatEnum pixel_format;
    ChannelTypeEnum channel_type;

    ImageFormat()
      : cols(0), rows(0), planes(0),
        pixel_format(VW_PIXEL_UNKNOWN),
        channel_type(VW_CHANNEL_UNKNOWN)
    {}

  };

  template <class PixelT>
  class ImageView {
  public:
    typedef PixelT pixel_type;
    ImageView();
    ImageView(unsigned cols, unsigned rows, unsigned planes=1);
    %extend {
      int get_cols() const { return self->cols(); }
      int get_rows() const { return self->rows(); }
      int get_planes() const { return self->planes(); }
      int get_channels() const { return self->channels(); }
      PixelT const& get_pixel( int x, int y, int p=0 ) { return self->operator()(x,y,p); }
      void set_pixel( PixelT const& val, int x, int y, int p=0 ) { self->operator()(x,y,p) = val; }
      ImageView get_region( int x, int y, int cols, int rows ) { return crop(*self,x,y,cols,rows); }
      void set_region( ImageView const& val, int x, int y, int cols, int rows ) { crop(*self,x,y,cols,rows) = val; }
      void set_region( PixelT const& val, int x, int y, int cols, int rows ) { fill( crop(*self,x,y,cols,rows), val ); }
      ImageView get_col( int index ) { return select_col(*self,index); }
      void set_col( ImageView const& val, int index ) { select_col(*self,index) = val; }
      void set_col( PixelT const& val, int index ) { fill( select_col(*self,index), val ); }
      ImageView get_row( int index ) { return select_row(*self,index); }
      void set_row( PixelT const& val, int index ) { fill( select_row(*self,index), val ); }
      void set_row( ImageView const& val, int index ) { select_row(*self,index) = val; }
      ImageView get_plane( int index ) { return select_plane(*self,index); }
      void set_plane( ImageView const& val, int index ) { select_plane(*self,index) = val; }
      void set_plane( PixelT const& val, int index ) { fill( select_plane(*self,index), val ); }
    }
    %pythoncode {
      __getitem__ = Image_getitem
      __setitem__ = Image_setitem
      cols = property(get_cols)
      rows = property(get_rows)
      planes = property(get_planes)
      channels = property(get_channels)
    }
  };

  template <class PixelT>
  class ImageViewRef {
  public:
    typedef PixelT pixel_type;
    ImageViewRef( ImageView<pixel_type> const& image );
    %extend {
      ImageViewRef() { return 0; }
      int get_cols() const { return self->cols(); }
      int get_rows() const { return self->rows(); }
      int get_planes() const { return self->planes(); }
      int get_channels() const { return self->channels(); }
      ImageView<PixelT> rasterize() const { return vw::copy(*self); }
      pixel_type get_pixel( int x, int y, int p=0 ) { return self->operator()(x,y,p); }
      ImageView<pixel_type> get_region( int x, int y, int cols, int rows ) { return crop(*self,x,y,cols,rows); }
      ImageView<pixel_type> get_col( int index ) { return select_col(*self,index); }
      ImageView<pixel_type> get_row( int index ) { return select_row(*self,index); }
      ImageView<pixel_type> get_plane( int index ) { return select_plane(*self,index); }
    }
    %pythoncode {
      __getitem__ = Image_getitem
      cols = property(get_cols)
      rows = property(get_rows)
      planes = property(get_planes)
      channels = property(get_channels)
      def ref(self):
        return self
    }
  };
}

%pythoncode {
  _pixel_image_table = dict()

  class Image(object):
    '''The standard Vision Workbench image class.'''

    def __init__( self, cols=0, rows=0, planes=1, ptype=None, pformat=None, ctype=None ):
      if isimage(cols):
        image = cols
        ptype = pixel._compute_pixel_type(image.pixel_type,ptype,pformat,ctype)
        if ptype != image.pixel_type:
          # we could cast here, but the easy way would create a cyclic dependency
          raise NotImplementedError, 'Incompatible pixel type'
        self.__dict__['impl'] = image.ref().rasterize()
      else:
        ptype = pixel._compute_pixel_type(pixel.PixelRGB_float32,ptype,pformat,ctype)
        self.__dict__['impl'] = _pixel_image_table[ptype](cols,rows,planes)

    def __getattr__(self,name):
      return getattr(self.impl,name)

    def __setattr__(self,name,value):
      return setattr(self.impl,name,value)

    __getitem__ = Image_getitem
    __setitem__ = Image_setitem
}

%define %instantiate_image_types(cname,ctype,pname,ptype,...)
  %template(ImageView_##pname) vw::ImageView<ptype >;
  %template(ImageViewRef_##pname) vw::ImageViewRef<ptype >;
  %pythoncode {
    _pixel_image_table[pixel.pname] = ImageView_##pname
    ImageView_##pname.pixel_type = pixel.pname
    ImageViewRef_##pname.pixel_type = pixel.pname
    try:
      ImageView_##pname.pixel_format = pixel.pname.pixel_format
      ImageViewRef_##pname.pixel_format = pixel.pname.pixel_format
    except:
      ImageView_##pname.pixel_format = pixel.PixelScalar
      ImageViewRef_##pname.pixel_format = pixel.PixelScalar
    ImageView_##pname.channel_type = pixel.cname
    ImageViewRef_##pname.channel_type = pixel.cname
    def _ImageView_##pname##_ref(self):
      return ImageViewRef_##pname(self)
    ImageView_##pname.ref = _ImageView_##pname##_ref
  }
%enddef

%instantiate_for_pixel_types(instantiate_image_types)


// Edge Extension

namespace vw {
  class ZeroEdgeExtension {};
  class ConstantEdgeExtension {};
  class PeriodicEdgeExtension {};
  class ReflectEdgeExtension {};
}

%inline %{
  template <class PixelT, class EdgeT>
  vw::ImageViewRef<PixelT> _edge_extend( vw::ImageViewRef<PixelT> const& image, long x_offset, long y_offset, long cols, long rows, EdgeT const& edge ) {
    return edge_extend( image, x_offset, y_offset, cols, rows, edge );
  }
%}

%define %instantiate_for_pixel_types_and_edge_extensions(macro,args...)
  %instantiate_for_pixel_types(macro,vw::ZeroEdgeExtension,args)
  %instantiate_for_pixel_types(macro,vw::ConstantEdgeExtension,args)
  %instantiate_for_pixel_types(macro,vw::PeriodicEdgeExtension,args)
  %instantiate_for_pixel_types(macro,vw::ReflectEdgeExtension,args)
%enddef

%define %instantiate_edge_extend(cname,ctype,pname,ptype,ee,...)
  %template(_edge_extend) _edge_extend<ptype, ee>;
%enddef

%instantiate_for_pixel_types_and_edge_extensions(instantiate_edge_extend)

%pythoncode {
  def edge_extend(image,xoffset=0,yoffset=0,cols=None,rows=None,edge=None):
    ref = image.ref()
    if cols is None: cols = ref.cols
    if rows is None: rows = ref.rows
    if edge is None: edge = ZeroEdgeExtension()
    return _edge_extend(ref,xoffset,yoffset,cols,rows,edge)
}
