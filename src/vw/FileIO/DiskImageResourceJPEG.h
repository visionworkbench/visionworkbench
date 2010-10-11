// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourceJPEG.h
///
/// Provides support for file formats via libJPEG.
///
#ifndef __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__
#define __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__

#include <cstdio>
#include <string>
#include <boost/shared_ptr.hpp>

#include <vw/FileIO/DiskImageResource.h>


namespace vw {

  class DiskImageResourceJPEG : public DiskImageResource {
  public:

    DiskImageResourceJPEG( std::string const& filename,
                           int subsample_factor = 1,
                           size_t byte_offset = 0)
      : DiskImageResource( filename )
    {
      m_quality = default_quality;
      m_file_ptr = NULL;
      m_jpg_compress_header = NULL;
      m_jpg_decompress_header = NULL;
      open( filename, subsample_factor, byte_offset );
    }

    DiskImageResourceJPEG( std::string const& filename,
                           ImageFormat const& format )
      : DiskImageResource( filename )
    {
      m_subsample_factor = default_subsampling_factor;
      m_quality = default_quality;
      m_file_ptr = NULL;
      m_jpg_compress_header = NULL;
      m_jpg_decompress_header = NULL;
      create( filename, format );
    }

    virtual ~DiskImageResourceJPEG();

    /// Returns the type of disk image resource.
    static std::string type_static() { return "JPEG"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    virtual void flush();

    /// Set the compression quality of the jpeg image.  The quality is
    /// a value between 0.0 and 1.0.  The lower the quality, the more
    /// lossy the compression.
    void set_quality(float quality) { m_quality = quality; }

    /// Set the default compression quality of jpeg images.
    static void set_default_quality(float quality) { default_quality = quality; }

    /// Set the subsample factor.  The default is no scaling.  Valid
    /// values are 1, 2, 4, and 8.  Smaller scaling ratios permit
    /// significantly faster decoding since fewer pixels need to be
    /// processed and a simpler IDCT method can be used.
    void set_subsample_factor(int subsample_factor) {
      // Cloes and reopen the file with the new subsampling factor
      flush();
      open(m_filename, subsample_factor);
    }

    /// Set the default subsampling factor
    static void set_default_subsample_factor(int subsample_factor) {
      default_subsampling_factor = subsample_factor;
    }

    void open( std::string const& filename,
               int subsample_factor = 1,
               size_t byte_offset = 0);

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

  private:
    // Forward declare an abstraction class that contains jpeg stuff.
    class vw_jpeg_decompress_context;
    friend class vw_jpeg_decompress_context;

    std::string m_filename;
    float m_quality;
    int m_subsample_factor;
    void* m_jpg_decompress_header;
    void* m_jpg_compress_header;
    FILE* m_file_ptr;
    size_t m_byte_offset;

    static int default_subsampling_factor;
    static float default_quality;

    /* The decompression context. We can't use an auto_ptr here because
     * vw_jpeg_decompress_context is forward-declared and the compiler
     * will not see its destructor, which means its destructor won't
     * actually get called :-( (See the compiler warnings when you try to
     * use std::auto_ptr<> here.) Boost's shared_ptr class, however, does
     * _not_ have this problem, and is perfectly safe to use in this case.
    */
    mutable boost::shared_ptr<vw_jpeg_decompress_context> ctx;

    /* Resets the decompression context and current point in the file to
     * the beginning.
    */
    void read_reset() const;

  };

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__
