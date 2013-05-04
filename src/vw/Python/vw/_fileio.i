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


%module fileio
%include "std_string.i"
%import "_image.i"

%{
#include <vw/Image.h>
#include <vw/FileIO.h>
%}

namespace vw {
  class DiskImageResource {
  public:
    virtual ~DiskImageResource();
    %extend {
      int get_cols() const { return self->cols(); }
      int get_rows() const { return self->rows(); }
      int get_planes() const { return self->planes(); }
      int get_channels() const { return self->channels(); }
      DiskImageResource( std::string const& filename ) { return vw::DiskImageResource::open( filename ); }
      PixelFormatEnum _pixel_format() const { return self->pixel_format(); }
      ChannelTypeEnum _channel_type() const { return self->channel_type(); }
    }
    %pythoncode {
      cols = property(get_cols)
      rows = property(get_rows)
      planes = property(get_planes)
      channels = property(get_channels)
      def get_pixel_format(self):
        return pixel._pixel_format_table[ self._pixel_format() ]
      def get_channel_type(self):
        return pixel._channel_type_table[ self._channel_type() ]
      def get_pixel_type(self):
        return self.pixel_format[self.channel_type]
      pixel_format = property(get_pixel_format)
      channel_type = property(get_channel_type)
      pixel_type = property(get_pixel_type)
    }
  };

  class DiskImageResourceJPEG : public DiskImageResource {
  public:
    virtual ~DiskImageResource();
    DiskImageResourceJPEG( std::string const& filename );
    static void set_default_quality( float quality );
  };

  class DiskImageResourcePNG : public DiskImageResource {
  public:
    DiskImageResourcePNG( std::string const& filename );
    unsigned num_comments() const;
    std::string const& get_comment_key  ( unsigned i ) const;
    std::string const& get_comment_value( unsigned i ) const;
    static void set_default_compression_level(int level);
  };
}

%inline %{
  template <class PixelT> void _read_image( vw::ImageView<PixelT>& image, vw::DiskImageResource& resource ) {
    vw::read_image( image, resource );
  }
  template <class PixelT> void _write_image( std::string const& filename, vw::ImageViewRef<PixelT> const& image ) {
    vw::write_image( filename, image );
  }
%}

namespace vw {
  template <class PixelT>
  class DiskImageView {
  public:
    typedef PixelT pixel_type;

    DiskImageView( std::string const& filename );
    DiskImageView( DiskImageResource *resource );

    virtual ~DiskImageView();

    %extend {
      int get_cols() const { return self->cols(); }
      int get_rows() const { return self->rows(); }
      int get_planes() const { return self->planes(); }
      int get_channels() const { return self->channels(); }
      std::string get_filename() const { return self->filename(); }
      vw::ImageViewRef<PixelT> ref() const { return *self; }
    }

    %pythoncode {
      cols = property(get_cols)
      rows = property(get_rows)
      planes = property(get_planes)
      channels = property(get_channels)
      filename = property(get_filename)
    }
  };
}

%pythoncode {
  def read_image(filename,ptype=None,pformat=None,ctype=None):
    resource = DiskImageResource(filename)
    im = image.Image(ptype=pixel._compute_pixel_type(resource.pixel_type(),ptype,pformat,ctype))
    _read_image(im.impl,resource)
    return im


  def write_image(filename,image):
    _write_image(filename,image.ref())


  class DiskImageView(object):
    '''A read-only view of an image on disk.'''

    _pixel_type_table = dict()

    def __init__(self,resource,ptype=None,pformat=None,ctype=None):
      if resource.__class__ is str:
        resource = DiskImageResource(resource)
      self.__dict__['pixel_type'] = pixel._compute_pixel_type(resource.pixel_type,ptype,pformat,ctype)
      resource.thisown = 0
      self.__dict__['_delegate'] = DiskImageView._pixel_type_table[self.pixel_type](resource)

    def __getattr__(self,name):
      return getattr(self._delegate,name)

    def __setattr__(self,name,value):
      return setattr(self._delegate,name,value)
}

%define %instantiate_fileio(cname,ctype,pname,ptype,...)
  %template(_read_image) _read_image<ptype >;
  %template(_write_image) _write_image<ptype >;
  %template(DiskImageView_##pname) vw::DiskImageView<ptype >;
  %pythoncode {
    DiskImageView._pixel_type_table[pixel.pname] = DiskImageView_##pname
  }
%enddef

%instantiate_for_pixel_types(instantiate_fileio)
