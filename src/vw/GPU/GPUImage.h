// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GPU_GPUIMAGE_H__
#define __VW_GPU_GPUIMAGE_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>

#include <vw/GPU/TexAlloc.h>
#include <vw/GPU/TexObj.h>
#include <vw/GPU/Expressions.h>
#include <vw/GPU/Interpolation.h>
#include <vw/GPU/EdgeExtension.h>


namespace vw {
namespace GPU {

  // Type is a dummy type representing 16-bit floating point precision
  // on a GPU, and is only for use as a template parameter to
  // GPUImage.
  class float16 {};

  //############################################################################
  //                          Class: GPUImageBase
  //############################################################################
  class GPUImageBase {
  protected:
    TexObj* m_texObj;
    int m_width;
    int m_height;
    int m_xOffset;
    int m_yOffset;
    bool m_isHomography;
    Matrix<float> m_homography;
    std::string m_interpolation_string;
    int m_interpolation_quality;
    EdgeExtensionType m_edge_extension_type;
  public:
    // Instance Functions - Ctors/Dtor
    GPUImageBase() {
      m_width = 0;
      m_height  = 0;
      m_texObj = NULL;
      m_isHomography = false;
    }

    GPUImageBase(int w, int h, Tex_Format format, Tex_Type type) {
      m_xOffset = 0;
      m_yOffset = 0;
      if(w ==0 || h == 0) {
        m_width = 0;
        m_height = 0;
        m_texObj = NULL;
        return;
      }
      m_texObj = TexAlloc::alloc(w, h, format, type);
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      m_isHomography = false;
    }

    GPUImageBase(int w, int h, Tex_Format format, Tex_Type type,
                 Tex_Format inputFormat, Tex_Type inputType, void* data)
    {
      m_xOffset = 0;
      m_yOffset = 0;
      if(w ==0 || h == 0) {
        m_width = 0;
        m_height = 0;
        m_texObj = NULL;
        return;
      }
      m_texObj = TexAlloc::alloc(w, h, format, type);
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      write(0, 0, m_width, m_height, inputFormat, inputType, data);
      m_isHomography = false;
    }

    ~GPUImageBase() {
      if(m_texObj) m_texObj->release();
    }


    GPUImageBase(const GPUImageBase& cpyTex) {
      m_width = cpyTex.m_width;
      m_height  = cpyTex.m_height;
      m_xOffset = cpyTex.m_xOffset;
      m_yOffset = cpyTex.m_yOffset;
      m_texObj = cpyTex.m_texObj;
      m_texObj->retain();
      m_isHomography = cpyTex.m_isHomography;
      m_homography = cpyTex.m_homography;
      m_interpolation_string = cpyTex.m_interpolation_string;
      m_interpolation_quality = cpyTex.m_interpolation_quality;
      m_edge_extension_type = cpyTex.m_edge_extension_type;

    }

    void copy_attributes(const GPUImageBase& cpyTex) {
      m_xOffset = 0;
      m_yOffset = 0;
      m_texObj = TexAlloc::alloc(cpyTex.width(), cpyTex.height(), cpyTex.format(), cpyTex.type());
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      m_isHomography = false;
    }

    // Instance Functions - Operators
    GPUImageBase& operator=(const GPUImageBase& cpyTex) {
      if(m_texObj)
        m_texObj->release();
      m_width = cpyTex.m_width;
      m_height  = cpyTex.m_height;
      m_xOffset = cpyTex.m_xOffset;
      m_yOffset = cpyTex.m_yOffset;
      m_texObj = cpyTex.m_texObj;
      m_texObj->retain();
      m_isHomography = cpyTex.m_isHomography;
      m_homography = cpyTex.m_homography;
      m_interpolation_string = cpyTex.m_interpolation_string;
      m_interpolation_quality = cpyTex.m_interpolation_quality;
      m_edge_extension_type = cpyTex.m_edge_extension_type;

      return *this;
    }
    // Instance Functions - ImageView Interface

    void set_size(int w, int h) { m_width = w; m_height = h; } // TEMP - Make Deep Operation

    int cols() const { return m_width; }

    int rows() const { return m_height; }

    void reset() {
      if(m_texObj)
        m_texObj->release();
      m_texObj = NULL;
      m_width = 0;
      m_height = 0;
    }

    int planes() { return 1; }

    // Instance Functions - Exclusive Interface

    template <class InterpT, class EdgeT>
    void apply_homography(Matrix<float>& homography, InterpT intepolation, EdgeT edge, unsigned int new_width = 0, unsigned int new_height = 0) {
      if(homography.cols() != 3 || homography.rows() != 3) {
        throw(Exception("[vw::GPU::GPUImage<PixelT>::apply_homography] Input Error: Homography Matrix must be 3x3"));
      }
       if(m_isHomography) {
        m_homography = homography * m_homography;
        int interpolation_quality = TraitsForInterpT<InterpT>::quality;
        if(interpolation_quality > m_interpolation_quality) {
          m_interpolation_quality = interpolation_quality;
          m_interpolation_string = TraitsForInterpT<InterpT>::ShaderString();
        }
      }
      else {
        m_homography.set_size(3, 3);
        m_homography.set_identity();
        m_homography = homography * m_homography;
        m_isHomography = true;
        m_interpolation_quality = TraitsForInterpT<InterpT>::quality;
        m_interpolation_string = TraitsForInterpT<InterpT>::ShaderString();
      }
      if(new_width) m_width = new_width;
      if(new_height) m_height = new_height;
      m_edge_extension_type = TraitsForEdgeT<EdgeT>::type;
    }

    void rasterize_homography() const;

    int width() const { return m_width; }
    int height() const { return m_height; }

    int real_width() const { return m_texObj->width(); }
    int real_height() const { return m_texObj->height(); }
    bool same_real_object(const GPUImageBase& other) { return false; }

    TexObj* GetTexObj() const { return m_texObj; } // TEMP - These 3 functions are only for use in the GPUImage copy constructor
    int OffsetX() const { return m_xOffset; }
    int OffsetY() const { return m_yOffset; }

    void translate(int x, int y) { m_xOffset += x; m_yOffset += y; }

    void write(Tex_Format inputFormat, Tex_Type inputType, void* data)
    { m_texObj->write(m_xOffset, m_yOffset, m_width, m_height, inputFormat, inputType, data); }

    void write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data)
    { m_texObj->write(x + m_xOffset, y + m_yOffset, w, h, inputFormat, inputType, data); }

    void read(Tex_Format outputFormat, Tex_Type outputType, void* data) const
    { rasterize_homography();m_texObj->read(m_xOffset, m_yOffset, m_width, m_height, outputFormat, outputType, data); }

    void read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data) const
    { rasterize_homography(); m_texObj->read(x + m_xOffset, y + m_yOffset, w, h, outputFormat, outputType, data); }

    GLuint target() const { return m_texObj->target(); }

    GLuint name() const { return m_texObj->name(); }

    void bind() const { m_texObj->bind(); }

    float x(float x) const { return m_texObj->x(x + m_xOffset); }

    float y(float y) const { return m_texObj->y(y + m_yOffset); }

    Tex_Format format() const { return m_texObj->format(); }

    Tex_Type type() const { return m_texObj->type(); }

    int num_channels() const {
      Tex_Format theFormat = format();
      if(theFormat == GPU_RGBA) return 4;
      else if(theFormat == GPU_RGB) return 3;
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
    static const Tex_Type gpu_type = GPU_FLOAT32;
    static const Tex_Type cpu_type = GPU_FLOAT32;
    typedef float input_type;
  };

  template <>
  struct TexTraitsForChannelT<float16> {
    static const Tex_Type gpu_type = GPU_FLOAT16;
    static const Tex_Type cpu_type = GPU_FLOAT32;
    typedef float input_type;
  };

  template <>
  struct TexTraitsForChannelT<vw::uint8> {
    static const Tex_Type gpu_type = GPU_UINT8;
    static const Tex_Type cpu_type = GPU_UINT8;
    typedef vw::uint8 input_type;
  };

  ///  Template Traits: TexTraitsForPixelT
  ///
  ///  (This abstraction is required in order to consider a
  ///  GPUImage<PixelRGBA<float16> > to be convertable to
  ///  ImageView<PixelRGBA<float> >. This is proper behavior since
  ///  float16 is not a real C++ type, but a stand in to represent the
  ///  16-bit precision floating point texture format that can be used
  ///  by the GPU.)
  template <class PixelT> struct TexTraitsForPixelT { typedef PixelT imageview_t; };

  template <> struct TexTraitsForPixelT<PixelRGBA<float16> > { typedef PixelRGBA<float> imageview_t; };
  template <> struct TexTraitsForPixelT<PixelRGB<float16> > { typedef PixelRGB<float> imageview_t; };
  template <> struct TexTraitsForPixelT<PixelGray<float16> > { typedef PixelGray<float> imageview_t; };
  template <> struct TexTraitsForPixelT<float16> { typedef float imageview_t; };

  //#################################################################################################################
  //                                 Proxy Class: GPUPixel<PixelT>
  //#################################################################################################################

  template <class PixelT> class GPUImage;

  template <class PixelT>
  class GPUPixel {
    TexObj* m_texObj;
    int m_xOffset;
    int m_yOffset;
  public:
    // Instance Functions - Ctors/Dtor
    GPUPixel() {
      m_texObj = NULL;
    }

    GPUPixel(const GPUImage<PixelT>& image, int xOffset, int yOffset) {
      m_xOffset = xOffset;
      m_yOffset = yOffset;
      m_texObj = image.m_texObj;
      m_texObj->retain();
    }

    ~GPUPixel() {
      if(m_texObj) m_texObj->release();
    }

    void operator=(const PixelT& pixel) {
      m_texObj->write(m_xOffset, m_yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(),
                     GPUImage<PixelT>::get_gpu_type_for_pixelt(), (void*) &pixel);
    }

    PixelT operator*() const {
      PixelT pixel;
      m_texObj->read(m_xOffset, m_yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(),
                    GPUImage<PixelT>::get_gpu_type_for_pixelt(), &pixel);
      return pixel;
    }


    operator PixelT() const {
      PixelT pixel;
      m_texObj->read(m_xOffset, m_yOffset, 1, 1, GPUImage<PixelT>::get_format_for_pixelt(),
                    GPUImage<PixelT>::get_gpu_type_for_pixelt(), &pixel);
      return pixel;
    }
  };

  //############################################################################
  //             Class: GPUImage<PixelT> (subclass of GPUImageBase)
  //############################################################################
  template <class PixelT>
  class GPUImage :  public GPUImageBase {
    friend class GPUPixel<PixelT>;
  public:

    static Tex_Format get_format_for_pixelt() {
      int num_channels = PixelNumChannels<PixelT>::value;
      Tex_Format tex_format;
      if(num_channels == 1)
        tex_format = GPU_RED;
      else if(num_channels == 4)
        tex_format = GPU_RGBA;
      else if(num_channels == 3)
        tex_format = GPU_RGB;
      else if(num_channels == 2) {
        printf("get_format_for_pixelt: 2 Channel not supported.\n");
        tex_format = GPU_RGB;
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
      m_width = 0;
      m_height = 0;
      m_xOffset = 0;
      m_yOffset = 0;
      m_texObj = NULL;
      m_isHomography = false;
    }

    GPUImage(int w, int h) {
      m_xOffset = 0;
      m_yOffset = 0;
      if(w ==0 || h == 0) {
        m_width = 0;
        m_height = 0;
        m_texObj = NULL;
        return;
      }
      Tex_Type gpu_type = get_gpu_type_for_pixelt();
      m_texObj = TexAlloc::alloc(w, h, get_format_for_pixelt(), gpu_type);
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      m_isHomography = false;
    }

    GPUImage(int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data) {
      m_width = w;
      m_height = h;
      m_xOffset = 0;
      m_yOffset = 0;
      m_isHomography = false;
      if(w ==0 || h == 0) {
        m_width = 0;
        m_height = 0;
        m_texObj = NULL;
        return;
      }
      Tex_Type gpu_type = get_gpu_type_for_pixelt();
      m_texObj = TexAlloc::alloc(w, h, get_format_for_pixelt(), gpu_type);
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      write(0, 0, m_width, m_height,  inputFormat, inputType, data);
    }

    GPUImage(const ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) {
      m_width = image.cols();
      m_height = image.rows();
      m_xOffset = 0;
      m_yOffset = 0;
      m_isHomography = false;
      Tex_Format tex_format = get_format_for_pixelt();
      m_texObj = TexAlloc::alloc(m_width, m_height, tex_format, get_gpu_type_for_pixelt());
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      write(0, 0, m_width, m_height, tex_format, get_cpu_type_for_pixelt(), &(image(0, 0)));
    }


    GPUImage(const GPUImage& cpyTex) {
      m_width = cpyTex.m_width;
      m_height  = cpyTex.m_height;
      m_isHomography = cpyTex.m_isHomography;
      m_homography = cpyTex.m_homography;
      m_interpolation_string = cpyTex.m_interpolation_string;
      m_interpolation_quality = cpyTex.m_interpolation_quality;
      m_edge_extension_type = cpyTex.m_edge_extension_type;
      m_xOffset = cpyTex.m_xOffset;
      m_yOffset = cpyTex.m_yOffset;
      m_texObj = cpyTex.m_texObj;
      m_texObj->retain();
    }

    GPUImage(const GPUImageBase& cpyTex) {
      if(cpyTex.format() != get_format_for_pixelt() || cpyTex.type() != get_gpu_type_for_pixelt()) {
        m_width = m_height = 0;
        throw(Exception("[vw::GPU::GPUImage<PixelT>::apply_homography] Input Error: Homography Matrix must be 3x3"));
      }
      *this = (GPUImage<PixelT>&) cpyTex;
    }

    GPUImage& operator=(const GPUImageBase& cpyTex) {
      if(cpyTex.format() != get_format_for_pixelt() || cpyTex.type() != get_gpu_type_for_pixelt()) {
        m_width = m_height = 0;
        throw(Exception("[vw::GPU::GPUImage<PixelT>::apply_homography] Input Error: Homography Matrix must be 3x3"));
      }
      *this = (GPUImage<PixelT>&) cpyTex;
      return *this;
    }

    GPUImage& operator=(const ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) {
      m_width = image.cols();
      m_height = image.rows();
      m_xOffset = 0;
      m_yOffset = 0;
      m_isHomography = false;
      Tex_Format tex_format = get_format_for_pixelt();
      m_texObj = TexAlloc::alloc(m_width, m_height, tex_format, get_gpu_type_for_pixelt());
      m_texObj->retain();
      m_width = m_texObj->width();
      m_height = m_texObj->height();
      write(0, 0, m_width, m_height, tex_format, get_cpu_type_for_pixelt(), &(image(0, 0)));
      return *this;
    }


    GPUImage& operator=(const GPUImage& cpyTex) {
      if(m_texObj)
        m_texObj->release();
      m_width = cpyTex.width();
      m_height  = cpyTex.height();
      m_isHomography = cpyTex.m_isHomography;
      m_homography = cpyTex.m_homography;
      m_interpolation_string = cpyTex.m_interpolation_string;
      m_interpolation_quality = cpyTex.m_interpolation_quality;
      m_edge_extension_type = cpyTex.m_edge_extension_type;
      m_xOffset = cpyTex.OffsetX();
      m_yOffset = cpyTex.OffsetY();
      m_texObj = cpyTex.m_texObj;
      m_texObj->retain();
      return *this;
    }

    const PixelT operator()(int col, int row) const {
      PixelT pixel;
      Tex_Type gpu_type = get_gpu_type_for_pixelt();
      Tex_Format tex_format = get_format_for_pixelt();
      read(col, row, 1, 1, tex_format, gpu_type, (void*) &pixel);
      return pixel;
    }

    GPUPixel<PixelT> pixel(int col, int row) {
      return GPUPixel<PixelT>(*this, col, row);
    }

    GPUPixel<PixelT> operator()(int col, int row) {
      return GPUPixel<PixelT>(*this, col, row);
    }

    void write_image_view(ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t>& image) const {
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
  inline void write_image(const std::string& filename, const GPUImage<PixelT>& image) {
    ImageView<typename TexTraitsForPixelT<PixelT>::imageview_t> imageView;
    image.write_image_view(imageView);
    write_image(filename, imageView);
  }



}} // namespace vw::GPU

#endif // __VW_GPU_GPUIMAGE_H__
