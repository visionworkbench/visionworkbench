%module fileio
%include "std_string.i"
%import "_image.i"

%{
#include <vw/Image.h>
#include <vw/FileiO.h>
%}

namespace vw {
  class DiskImageResource {
  public:
    virtual ~DiskImageResource();
    int32 cols() const;
    int32 rows() const;
    int32 planes() const;
    int32 channels() const;
    %extend {
      DiskImageResource( std::string const& filename ) { return vw::DiskImageResource::open( filename ); }
      PixelFormatEnum _pixel_format() const { return self->pixel_format(); }
      ChannelTypeEnum _channel_type() const { return self->channel_type(); }
    }
    %pythoncode {
      def pixel_format(self):
        return pixel._pixel_format_table[ self._pixel_format() ]
      def channel_type(self):
        return pixel._channel_type_table[ self._channel_type() ]
    }
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

    virtual ~DiskImageView();
    
    unsigned cols() const;
    unsigned rows() const;
    unsigned planes() const;
    unsigned channels() const;

    std::string filename() const;

    %extend {
      vw::ImageViewRef<PixelT> ref() const { return *self; }
    }
  };
}

%pythoncode {
  def _compute_pixel_type(r,ptype,pformat,ctype):
    if ptype is None:
      if ctype is None: ctype = r.channel_type()
      if pformat is None: pformat = r.pixel_format()
      ptype = pformat[ctype]
    elif pformat is not None or ctype is not None:
      raise Exception, "Cannot specify both ptype and pformat/ctype"
    return ptype

	    
  def read_image(filename,ptype=None,pformat=None,ctype=None):
    r = DiskImageResource(filename)
    im = image.Image(ptype=_compute_pixel_type(r,ptype,pformat,ctype))
    _read_image( im, r )
    return im


  def write_image(filename,image):
    _write_image(filename,image.ref())


  class DiskImageView(object):
    '''A read-only view of an image on disk.'''

    _pixel_type_table = dict()

    def __init__(self,filename,ptype=None,pformat=None,ctype=None):
      r = DiskImageResource(filename)
      self.pixel_type = _compute_pixel_type(r,ptype,pformat,ctype)
      self.impl = self._pixel_type_table[self.pixel_type](filename)
    
    def _get_cols(self):
      return self.impl.cols()
    cols = property(_get_cols)

    def _get_rows(self):
      return self.impl.rows()
    rows = property(_get_rows)

    def _get_planes(self):
      return self.impl.planes()
    planes = property(_get_planes)

    def _get_channels(self):
      return self.impl.channels()
    channels = property(_get_channels)

    def _get_filename(self):
      return self.impl.filename()
    filename = property(_get_filename)

    def ref(self):
      return self.impl.ref() 
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
