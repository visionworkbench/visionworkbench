
#ifndef TEXREF_H
#define TEXREF_H

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>


#include <vw/GPU/TexAlloc.h>
#include <vw/GPU/TexBlock.h>
#include <vw/GPU/TexObj.h>
#include <vw/GPU/Expressions.h>
#include <vw/GPU/Interpolation.h>
#include <vw/GPU/EdgeExtension.h>


namespace vw { namespace GPU {

  class float16 {  // Type is a dummy type, and is only for use as a template parameter to GPUImage.
  private:
    float16() { }
  }; 

//#################################################################################################################
//                                               Class: GPUImageBase
//#################################################################################################################

class GPUImageBase : public ShaderNode_Base {
 protected:
// Member Variables - Protected
  TexBlock* _texBlock;
  int _width;
  int _height;
  int _xOffset;
  int _yOffset;
  bool _isHomography;
  Matrix<float> _homography;
  string _interpolation_string;
  int _interpolation_quality;
  EdgeExtensionType _edge_extension_type;
  public:
// Instance Functions - Ctors/Dtor
  GPUImageBase() { 
    _width = 0;
    _height  = 0;
    _texBlock = NULL;
    _isHomography = false;
  }

  GPUImageBase(int w, int h, Tex_Format format, Tex_Type type) {
    _xOffset = 0;
    _yOffset = 0;
    if(w ==0 || h == 0) {
	  _width = 0;
	  _height = 0;
	  _texBlock = NULL;
      return;
	}
    _texBlock = TexAlloc::alloc(w, h, format, type);
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    _isHomography = false;
  }

  GPUImageBase(int w, int h, Tex_Format format, Tex_Type type, 
		 Tex_Format inputFormat, Tex_Type inputType, void* data) 
  {
    _xOffset = 0; 
    _yOffset = 0;
    if(w ==0 || h == 0) {
	  _width = 0;
	  _height = 0;
	  _texBlock = NULL;
      return;
	}
    _texBlock = TexAlloc::alloc(w, h, format, type);
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    write(0, 0, _width, _height, inputFormat, inputType, data);
    _isHomography = false; 
  }

  ~GPUImageBase() {
    if(_texBlock) _texBlock->release(); 
  }

  GPUImageBase(const GPUImageBase& cpyTex) {
    _width = cpyTex._width;
    _height  = cpyTex._height;
    _xOffset = cpyTex._xOffset;
    _yOffset = cpyTex._yOffset;
    _texBlock = cpyTex._texBlock;
    _texBlock->retain();
    _isHomography = cpyTex._isHomography;
    _homography = cpyTex._homography;
    _interpolation_string = cpyTex._interpolation_string;
    _interpolation_quality = cpyTex._interpolation_quality;
    _edge_extension_type = cpyTex._edge_extension_type;

  }
// Instance Functions - Operators
  GPUImageBase& operator=(const GPUImageBase& cpyTex) {
    if(_texBlock)
      _texBlock->release();
    _width = cpyTex._width;
    _height  = cpyTex._height;
    _xOffset = cpyTex._xOffset;
    _yOffset = cpyTex._yOffset;
    _texBlock = cpyTex._texBlock;
    _texBlock->retain();
    _isHomography = cpyTex._isHomography;
    _homography = cpyTex._homography;
    _interpolation_string = cpyTex._interpolation_string;
    _interpolation_quality = cpyTex._interpolation_quality;
    _edge_extension_type = cpyTex._edge_extension_type;

    return *this;
  }
// Instance Functions - ImageView Interface

  void set_size(int w, int h) { _width = w; _height = h; } // TEMP - Make Deep Operation

  int cols() const { return _width; }
  
  int rows() const { return _height; }

  void reset() { 
    if(_texBlock)
      _texBlock->release();
    _texBlock = NULL;
    _width = 0;
    _height = 0;
  }
  
  int planes() { return 1; }
  
// Instance Functions - Exclusive Interface

  template <class InterpT, class EdgeT>
  void apply_homography(Matrix<float>& homography, InterpT intepolation, EdgeT edge, unsigned int new_width = 0, unsigned int new_height = 0) {
    if(homography.cols() != 3 || homography.rows() != 3) {
      throw(Exception("[vw::GPU::GPUImage<PixelT>::apply_homography] Input Error: Homography Matrix must be 3x3"));
    }
    if(_isHomography) {
      _homography = homography * _homography;
      int interpolation_quality = TraitsForInterpT<InterpT>::quality;
      if(interpolation_quality > _interpolation_quality) { 
	_interpolation_quality = interpolation_quality;
	_interpolation_string = TraitsForInterpT<InterpT>::ShaderString();
      }
    }
    else {
      _homography.set_size(3, 3);
      _homography.set_identity();
      _homography = homography * _homography;      
      _isHomography = true;
      _interpolation_quality = TraitsForInterpT<InterpT>::quality;
      _interpolation_string = TraitsForInterpT<InterpT>::ShaderString();
    }
    if(new_width) _width = new_width;
    if(new_height) _height = new_height;
    _edge_extension_type = TraitsForEdgeT<EdgeT>::type;
  }

  void rasterize_homography() const;

  int width() const { return _width; }
  int height() const { return _height; }

  int real_width() const { return _texBlock->real_width(); }
  int real_height() const { return _texBlock->real_height(); }
  bool same_real_object(const GPUImageBase& other) { return _texBlock->same_real_object(*(other._texBlock)); }
 
   TexBlock* GetTexBlock() const { return _texBlock; } // TEMP - These 3 functions are only for use in the GPUImage copy constructor
  int OffsetX() const { return _xOffset; }
  int OffsetY() const { return _yOffset; }

  void translate(int x, int y) { _xOffset += x; _yOffset += y; }

  void write(Tex_Format inputFormat, Tex_Type inputType, void* data) 
    { _texBlock->write(_xOffset, _yOffset, _width, _height, inputFormat, inputType, data); }

  void write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data) 
    { _texBlock->write(x + _xOffset, y + _yOffset, w, h, inputFormat, inputType, data); }

  void read(Tex_Format outputFormat, Tex_Type outputType, void* data)
    { rasterize_homography();_texBlock->read(_xOffset, _yOffset, _width, _height, outputFormat, outputType, data); }

  void read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data)
    { rasterize_homography(); _texBlock->read(x + _xOffset, y + _yOffset, w, h, outputFormat, outputType, data); }

  GLuint target() const { return _texBlock->target(); }

  GLuint name() const { return _texBlock->name(); }

  void bind() const { _texBlock->bind(); }

  float x(float x) const { return _texBlock->x(x + _xOffset); }

  float y(float y) const { return _texBlock->y(y + _yOffset); }

  Tex_Format format() const { return _texBlock->format(); }

  Tex_Type type() const { return _texBlock->type(); }

  int num_channels() const {
    Tex_Format theFormat = format();
    if(theFormat == TEX_RGBA) return 4;
    else if(theFormat == TEX_RGB) return 3;
    else return 1;
  }
};


//#################################################################################################################
//                                      Template Traits: TexTraitsForChannelT
//#################################################################################################################

 template <class PixelT> 
 struct TexTraitsForChannelT { };

 template <> 
 struct TexTraitsForChannelT<float> {
   static const Tex_Type gpu_type = TEX_FLOAT32;
   static const Tex_Type cpu_type = TEX_FLOAT32;
   typedef float input_type;
 };

 template <> 
 struct TexTraitsForChannelT<float16> {
   static const Tex_Type gpu_type = TEX_FLOAT16;
   static const Tex_Type cpu_type = TEX_FLOAT32;
   typedef float input_type;
 };

 template <> 
 struct TexTraitsForChannelT<vw::uint8> {
   static const Tex_Type gpu_type = TEX_UINT8;
   static const Tex_Type cpu_type = TEX_UINT8;
   typedef vw::uint8 input_type;
 };

//#################################################################################################################
//                                      Template Traits: TexTraitsForPixelT
//   (This abstraction is required in order to consider a GPUImage<PixelRGBA<float16> > to be convertable to
//    ImageView<PixelRGBA<float> >. This is proper behavior since float16 is not a real C++ type, but a stand in to
//    represent the 16-bit precision floating point texture format that can be used by the GPU.)
//#################################################################################################################

template <class PixelT>
struct TexTraitsForPixelT { };

 
 template <> struct TexTraitsForPixelT<PixelRGBA<float> > { typedef PixelRGBA<float> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelRGB<float> > { typedef PixelRGB<float> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelGray<float> > { typedef PixelGray<float> imageview_t; };
 template <> struct TexTraitsForPixelT<float> { typedef float imageview_t; };

 template <> struct TexTraitsForPixelT<PixelRGBA<float16> > { typedef PixelRGBA<float> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelRGB<float16> > { typedef PixelRGB<float> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelGray<float16> > { typedef PixelGray<float> imageview_t; };
 template <> struct TexTraitsForPixelT<float16> { typedef float imageview_t; };

 template <> struct TexTraitsForPixelT<PixelRGBA<vw::uint8> > { typedef PixelRGBA<vw::uint8> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelRGB<vw::uint8> > { typedef PixelRGB<vw::uint8> imageview_t; };
 template <> struct TexTraitsForPixelT<PixelGray<vw::uint8> > { typedef PixelGray<vw::uint8> imageview_t; };
 template <> struct TexTraitsForPixelT<vw::uint8> { typedef vw::uint8 imageview_t; };

//#################################################################################################################
//                                 Proxy Class: GPUPixel<PixelT>
//#################################################################################################################

template <class PixelT> class GPUImage;

template <class PixelT>
class GPUPixel { 
  TexBlock* _texBlock;
  int _xOffset;
  int _yOffset;
  public:
// Instance Functions - Ctors/Dtor
  GPUPixel() { 
    _texBlock = NULL;
  }

  GPUPixel(const GPUImage<PixelT>& image, int xOffset, int yOffset) {
    _xOffset = xOffset;
    _yOffset = yOffset;
    _texBlock = image._texBlock;
    _texBlock->retain();
  }

  ~GPUPixel() {
    if(_texBlock) _texBlock->release(); 
  }
  
  void operator=(const PixelT& pixel) {
	_texBlock->write(_xOffset, _yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(), 
					  GPUImage<PixelT>::get_gpu_type_for_pixelt(), (void*) &pixel);
  }
  
  PixelT operator*() {
	PixelT pixel;
	_texBlock->read(_xOffset, _yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(), 
					  GPUImage<PixelT>::get_gpu_type_for_pixelt(), &pixel);
	return pixel;
  }
  
  
  operator PixelT() {
	PixelT pixel;
	_texBlock->read(_xOffset, _yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(), 
					  GPUImage<PixelT>::get_gpu_type_for_pixelt(), &pixel);
	return pixel;
  }

};

//#################################################################################################################
//                                 Class: GPUImage<PixelT> (subclass of GPUImageBase)
//#################################################################################################################

template <class PixelT>
class GPUImage :  public GPUImageBase { // public ImageViewBase<GPUImage<PixelT> >, {
 friend class GPUPixel<PixelT>;
 public:
  static Tex_Format get_format_for_pixelt() {
    int num_channels = PixelNumChannels<PixelT>::value;
    Tex_Format tex_format;
    if(num_channels == 1)
      tex_format = TEX_R;
    else if(num_channels == 4)
      tex_format = TEX_RGBA;
    else if(num_channels == 3)
      tex_format = TEX_RGB;
    else if(num_channels == 2) {
      printf("get_format_for_pixelt: 2 Channel not supported.\n");
      tex_format = TEX_RGB;
    }
    return tex_format;
  }

  static Tex_Type get_gpu_type_for_pixelt() {
    return TexTraitsForChannelT<typename PixelChannelType<PixelT>::type >::gpu_type;
  }

  static Tex_Type get_cpu_type_for_pixelt() {
    return TexTraitsForChannelT<typename PixelChannelType<PixelT>::type >::cpu_type;
  }

  GPUImage() {
    _width = 0;
    _height = 0;
    _xOffset = 0;
    _yOffset = 0;
    _texBlock = NULL;
    _isHomography = false;
  }
  
  GPUImage(int w, int h) {
    _xOffset = 0;
    _yOffset = 0;
    if(w ==0 || h == 0) {
	  _width = 0;
	  _height = 0;
	  _texBlock = NULL;
      return;
	}
    Tex_Type gpu_type = get_gpu_type_for_pixelt();
    _texBlock = TexAlloc::alloc(w, h, get_format_for_pixelt(), gpu_type);
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    _isHomography = false;
  }
   
  GPUImage(int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data) {
    _width = w;
    _height = h;
    _xOffset = 0; 
    _yOffset = 0;
    _isHomography = false;
    if(w ==0 || h == 0) {
	  _width = 0;
	  _height = 0;
	  _texBlock = NULL;
      return;
	}
    Tex_Type gpu_type = get_gpu_type_for_pixelt();
    _texBlock = TexAlloc::alloc(w, h, get_format_for_pixelt(), gpu_type);
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    write(0, 0, _width, _height,  inputFormat, inputType, data); 
  }

  GPUImage(const ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) {
    _width = image.cols();
    _height = image.rows();
    _xOffset = 0;
    _yOffset = 0;
    _isHomography = false;
    Tex_Format tex_format = get_format_for_pixelt();
    _texBlock = TexAlloc::alloc(_width, _height, tex_format, get_gpu_type_for_pixelt());
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    write(0, 0, _width, _height, tex_format, get_cpu_type_for_pixelt(), &(image(0, 0)));
  }  


  GPUImage(const GPUImage& cpyTex) {
    _width = cpyTex._width;
    _height  = cpyTex._height;
    _isHomography = cpyTex._isHomography;
    _homography = cpyTex._homography;
    _interpolation_string = cpyTex._interpolation_string;
    _interpolation_quality = cpyTex._interpolation_quality;
    _edge_extension_type = cpyTex._edge_extension_type;
    _xOffset = cpyTex._xOffset;
    _yOffset = cpyTex._yOffset;
    _texBlock = cpyTex._texBlock;
    _texBlock->retain();
  }

  GPUImage(const GPUImageBase& cpyTex) {
	if(cpyTex.format() != get_format_for_pixelt() || cpyTex.type() != get_gpu_type_for_pixelt()) {
		_width = _height = 0;
		throw(Exception("[vw::GPU::GPUImage<PixelT>::apply_homography] Input Error: Homography Matrix must be 3x3"));
    }
	*this = (GPUImage<PixelT>&) cpyTex;
    /*
    _width = cpyTex.width();
    _height  = cpyTex.height();
    _xOffset = cpyTex.OffsetX();
    _yOffset = cpyTex.OffsetY();
    _texBlock = cpyTex.GetTexBlock();
    _texBlock->retain();
    */
  }
  
  GPUImage& operator=(const ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) {
    _width = image.cols();
    _height = image.rows();
    _xOffset = 0;
    _yOffset = 0;
    _isHomography = false;
    Tex_Format tex_format = get_format_for_pixelt();
    _texBlock = TexAlloc::alloc(_width, _height, tex_format, get_gpu_type_for_pixelt());
    _texBlock->retain();
    _width = _texBlock->width();
    _height = _texBlock->height();
    write(0, 0, _width, _height, tex_format, get_cpu_type_for_pixelt(), &(image(0, 0)));
    return *this;
  }


  GPUImage& operator=(const GPUImage& cpyTex) {
    if(_texBlock)
      _texBlock->release();
    _width = cpyTex.width();
    _height  = cpyTex.height();
    _isHomography = cpyTex._isHomography;
    _homography = cpyTex._homography;
    _interpolation_string = cpyTex._interpolation_string;
    _interpolation_quality = cpyTex._interpolation_quality;
    _edge_extension_type = cpyTex._edge_extension_type;
    _xOffset = cpyTex.OffsetX();
    _yOffset = cpyTex.OffsetY();
    _texBlock = cpyTex._texBlock;
    _texBlock->retain();
    return *this;
  }

  const PixelT operator()(int col, int row) {
    PixelT pixel;
    Tex_Type gpu_type = get_gpu_type_for_pixelt();
    Tex_Format tex_format = get_format_for_pixelt();
    read(col, row, 1, 1, tex_format, gpu_type, (void*) &pixel);
    return pixel;
  }

  GPUPixel<PixelT> pixel(int col, int row) {
	return GPUPixel<PixelT>(*this, col, row);
  }

  void write_image_view(ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) {
    Tex_Type cpu_type =  get_cpu_type_for_pixelt(); 
    image.set_size(width(), height());
    read(get_format_for_pixelt(), cpu_type, &(image(0, 0)));
  }

  // void rasterize(ImageView<PixelT>& image) {
  // Tex_Type gpu_type =  get_gpu_type_for_pixelt(); 
  // read(get_format_for_pixelt(), gpu_type, &(image(0, 0)));
  // }

};

//#################################################################################################################
//                      Free Functions: read_image and write_image (overloads for vw/FileIO.h)
//#################################################################################################################

  template <class PixelT>
  inline void read_image(GPUImage<PixelT>& image, const std::string& filename) {
    ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t> imageView;
    read_image(imageView, filename);
    image = imageView;
  }

  template <class PixelT>
  inline void write_image(const std::string& filename, GPUImage<PixelT>& image) {
    ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t> imageView;
    image.write_image_view(imageView);
    write_image(filename, imageView);
  }



} } // namespaces GPU, vw

#endif
