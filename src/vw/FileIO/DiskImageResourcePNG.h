// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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


/// \file DiskImageResourcePNG.h
///
/// Provides support for the PNG file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCEPNG_H__
#define __VW_FILEIO_DISKIMAGERESOURCEPNG_H__

#include <vw/vw_config.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/PixelTypes.h>
#include <string>
#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {

  class DiskImageResourcePNG : public DiskImageResource {
  public:

    // Forward declare the options structure, which gets defined later on.
    struct Options;

    // The standard DiskImageResource interface:

    DiskImageResourcePNG( std::string const& filename );

    DiskImageResourcePNG( std::string const& filename,
                          ImageFormat const& format );

    // Additional constructor for passing in options.
    DiskImageResourcePNG( std::string const& filename,
                          ImageFormat const& format,
                          Options const& options );

    virtual ~DiskImageResourcePNG();

    /// Returns the type of disk image resource.
    static std::string type_static() { return "PNG"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read(ImageBuffer const& buf, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    void create( std::string const& filename,
                 ImageFormat const& format,
                 Options const& options );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    // The PNG-specific interface:

    struct Comment {
      std::string key, text, lang, lang_key;
      bool utf8, compressed;
    };

    unsigned num_comments() const;
    Comment const& get_comment( unsigned i ) const;
    std::string const& get_comment_key  ( unsigned i ) const;
    std::string const& get_comment_value( unsigned i ) const;
    void write_comments( std::vector<Comment> const& comments );

    // Options that can be passed in while writing. The default is
    // generally 'good enough'.
    struct Options {
      int compression_level; // From 0 to 9
      bool using_interlace; // Interlace output.
      bool using_palette; // Output as a palette.
      bool using_palette_indices; // For setting your own palette.
      bool using_palette_alpha; // If your palette has an alpha channel.
      std::vector<PixelRGBA<uint8>> palette; // The palette, set manually.

      Options();
    };

    static void set_default_compression_level(int level);

    // Convenience functions:

    template <typename ImageT>
    static void write_palette_file(std::string const& filename,
                                   ImageT const& image,
                                   std::vector<PixelRGBA<uint8>> const& palette) {
      Options options;
      options.using_palette = true;
      options.using_palette_indices = true;
      options.using_palette_alpha = true;
      options.palette = palette;
      DiskImageResourcePNG png(filename, image.format(), options);
      png.write(image.buffer(), BBox2i(0, 0, image.cols(), image.rows()));
    }

    virtual bool has_block_write()  const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read()   const {return true;}
    virtual bool has_nodata_read()  const {return false;}

    virtual Vector2i block_read_size() const { return m_block_size; }

  private:

    // vw_png_context is declared in the cc file, and unused elsewhere,
    // but we have a pointer to it. So forward declare.
    struct vw_png_context;
    struct vw_png_read_context;
    struct vw_png_write_context;

    // Block reading is supported. Block writing is not
    Vector2i m_block_size;

    mutable boost::shared_ptr<vw_png_context> m_ctx;

    void read_reset() const;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOURCEPNG_H__
