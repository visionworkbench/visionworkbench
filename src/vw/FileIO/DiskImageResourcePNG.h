// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file DiskImageResourcePNG.h
/// 
/// Provides support for the PNG file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEPNG_H__
#define __VW_FILEIO_DISKIMAGERESOUCEPNG_H__

#include <vw/config.h>

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>

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

    // Options that can be passed in while writing. The default is
    // generally 'good enough'.
    struct Options {
      int compression_level; // From 0 to 9
      bool using_interlace; // Interlace output.
      bool using_palette; // Output as a palette.
      bool using_palette_indices; // For setting your own palette.
      bool using_palette_alpha; // If your palette has an alpha channel.
      ImageView<PixelRGBA<uint8> > palette; // The palette, set manually.

      Options();
    };

    static void set_default_compression_level(int level);

    // Convenience functions:

    static void write_palette_file( std::string const& filename,
                                    ImageView<uint8> const& image,
                                    ImageView<PixelRGBA<uint8> > const& palette )
    {
      Options options;
      options.using_palette = true;
      options.using_palette_indices = true;
      options.using_palette_alpha = true;
      options.palette = palette;

      DiskImageResourcePNG png( filename, image.format(), options );
      png.write( image.buffer(), BBox2i(0,0,image.cols(),image.rows()) );
    }

  private:
    // vw_png_context is declared in the cc file, and unused elsewhere, 
    // but we have a pointer to it. So forward declare.
    struct vw_png_context;
    struct vw_png_read_context;
    struct vw_png_write_context;

    mutable boost::shared_ptr<vw_png_context> m_ctx;

    void read_reset() const;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEPNG_H__
