%module mosaic
%include "std_string.i"
%include "std_vector.i"
%import "_image.i"

%{
#include <vw/Image.h>
#include <vw/Mosaic.h>
%}

%inline %{
  struct PatchInfo {
    std::string filename;
    unsigned level;
    vw::BBox2i image_bbox;
    vw::BBox2i visible_bbox;
    vw::BBox2i region_bbox;
  };
%}

namespace vw {
namespace mosaic {

  template <class PixelT>
  class ImageQuadTreeGenerator {
  public:
    virtual ~ImageQuadTreeGenerator();
    typedef ::PatchInfo PatchInfo;

    ImageQuadTreeGenerator( std::string const& tree_name, vw::ImageViewRef<PixelT> const& source );
    
    void generate();

    void set_crop_bbox( vw::BBox2i const& bbox );
    vw::BBox2i const& get_crop_bbox() const;
    void set_output_image_file_type( std::string const& extension );
    std::string const& get_output_image_file_type() const;
    void set_patch_size( unsigned size );
    unsigned get_patch_size() const;
    void set_patch_overlap( unsigned overlap );
    unsigned get_patch_overlap() const;
    void set_levels_per_directory( unsigned levels_per_directory );
    unsigned get_levels_per_directory() const;
    void set_crop_images( bool crop );
    bool get_crop_images() const;
    void set_base_dir( std::string const& base_dir );
    std::string const& get_base_dir() const;
  };

} // namespace mosaic
} // namespace vw

%pythoncode {

  class ImageQuadTreeGenerator(object):
    '''The base quad-tree generator type.'''
    
    _pixel_type_table = dict()

    def __init__(self,name,source,**args):
      self.impl = self._pixel_type_table[source.pixel_type](name,source.ref())
      for (key,value) in args.iteritems():
        self.__class__.__dict__[key].fset(self,value)

    def generate(self):
      self.impl.generate()

    def _get_crop_bbox(self):
      return self.impl.get_crop_bbox()
    def _set_crop_bbox(self,bbox):
      self.impl.set_crop_bbox(bbox)
    crop_bbox = property(_get_crop_bbox,_set_crop_bbox)

    def _get_output_image_file_type(self):
      return self.impl.get_output_image_file_type()
    def _set_output_image_file_type(self,file_type):
      self.impl.set_output_image_file_type(file_type)
    output_image_file_type = property(_get_output_image_file_type,_set_output_image_file_type)

    def _get_patch_size(self):
      return self.impl.get_patch_size()
    def _set_patch_size(self,size):
      self.impl.set_patch_size(size)
    patch_size = property(_get_patch_size,_set_patch_size)

    def _get_patch_overlap(self):
      return self.impl.get_patch_overlap()
    def _set_patch_overlap(self,overlap):
      self.impl.set_patch_overlap(overlap)
    patch_overlap = property(_get_patch_overlap,_set_patch_overlap)

    def _get_levels_per_directory(self):
      return self.impl.get_levels_per_directory()
    def _set_levels_per_directory(self,levels):
      self.impl.set_levels_per_directory(levels)
    levels_per_directory = property(_get_levels_per_directory,_set_levels_per_directory)
    
    def _get_crop_images(self):
      return self.impl.get_crop_images()
    def _set_crop_images(self,crop):
      self.impl.set_crop_images(crop)
    crop_images = property(_get_crop_images,_set_crop_images)

    def _get_base_dir(self):
      return self.impl.get_base_dir()
    def _set_base_dir(self,dir):
      self.impl.set_base_dir(dir)
    base_dir = property(_get_base_dir,_set_base_dir)
}

%define %instantiate_mosaic_types(cname,ctype,pname,ptype,...)
  %template(ImageQuadTreeGenerator_##pname) vw::mosaic::ImageQuadTreeGenerator<ptype >;
  %pythoncode {
    ImageQuadTreeGenerator._pixel_type_table[pixel.pname] = ImageQuadTreeGenerator_##pname
  }
%enddef
   
%instantiate_for_pixel_types(instantiate_mosaic_types)
