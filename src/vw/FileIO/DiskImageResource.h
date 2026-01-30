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


/// \file FileIO/DiskImageResource.h
///
/// An abstract base class referring to an image on disk.
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCE_H__
#define __VW_FILEIO_DISKIMAGERESOURCE_H__

#include <set>
#include <string>

// Can't do much about warnings in boost except to hide them
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/type_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <vw/Core/Features.h>
#include <vw/Core/Log.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageIO.h>

namespace vw {

  // Return a smart pointer, this is easier to manage
  class DiskImageResource;
  boost::shared_ptr<DiskImageResource> DiskImageResourcePtr(std::string const& image_file);
  
  // *******************************************************************
  // The DiskImageResource abstract base class
  // *******************************************************************

  /// A base class for ImageResource objects where the buffer is on disk.
  /// - Base class from which specific file handlers derive.
  /// - Noncopyable because every impl is noncopyable
  class DiskImageResource : public ImageResource,
                            private boost::noncopyable {
  public:

    virtual ~DiskImageResource() {};

    /// Returns the type of disk image resource. Subclasses should
    /// implement this by calling a static function type_static().
    virtual std::string type() = 0;

    /// Returns the channel type of an image on disk.
    virtual ImageFormat format() const { return m_format; }

    /// Return the filename of the disk image file.
    std::string filename() const { return m_filename; }

    /// Create a new DiskImageResource of the appropriate type
    /// pointing to an existing file on disk.
    ///
    /// Don't forget to delete the DiskImageResource object when
    /// you're finished with it!
    static DiskImageResource* open( std::string const& filename );

    /// Create a new DiskImageResource of the appropriate type
    /// pointing to a newly-created empty file on disk.
    ///
    /// The underlying driver chooses a file format that it supports
    /// that matches the one you requested as closely as possible.
    /// If you care exactly what format it chose, you can ask the
    /// resource after it's created.  If you need finer-grained
    /// control you must manually create a resource of the appropraite
    /// type.  Don't forget to delete this DiskImageResource object
    /// when you're finished with it!
    static DiskImageResource* create( std::string const& filename, ImageFormat const& format );
    static DiskImageResource* create( std::string const& filename, ImageFormat const& format, std::string const& type );

    typedef DiskImageResource* (*construct_open_func)( std::string const& filename );

    typedef DiskImageResource* (*construct_create_func)( std::string const& filename,
                                                         ImageFormat const& format );

    /// Specify that a certain file extension should be opened with a certain
    ///  DiskImageResource derived type.
    static void register_file_type( std::string const& extension,
                                    std::string const& disk_image_resource_type,
                                    construct_open_func   open_func,
                                    construct_create_func create_func );

    /// This method is called automatically the first time you read or
    /// write a file, however in some situations you may want to
    /// register the default file types _before_ registering your own
    /// file types so that you can override the default FileIO drivers.

    // Specify whether values should be rescaled when converting
    // from one channel type to another during reads or writes.
    // Defaults to the value of default_rescale, which defaults to true.
    void set_rescale(bool rescale);

    // Specify a global default for rescale.  This is a dangerous
    // function to use, because some higher-level Vision Workbench
    // features rely on rescaling.  Use at your own risk.
    static void set_default_rescale(bool rescale);

    // TODO: This has always been the default, but it probably shouldn't be.
    virtual void flush() {}

  protected:
    DiskImageResource( std::string const& filename ) : m_filename(filename), m_rescale(default_rescale) {}
    ImageFormat m_format;
    std::string m_filename;
    bool        m_rescale;
    static bool default_rescale;
  };


  // *******************************************************************
  // Free functions using the DiskImageResource interface
  // *******************************************************************

  /// Read an image on disk into a vw::ImageView<T> object.
  template <class PixelT>
  void read_image( ImageView<PixelT>& in_image, const std::string &filename ) {

    VW_OUT(InfoMessage, "fileio") << "\tLoading image: " << filename << "\t";

    // Open the file for reading
    boost::shared_ptr<DiskImageResource> r(DiskImageResourcePtr(filename));

    VW_OUT(InfoMessage, "fileio") << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

    // Read the data
    read_image(in_image, *r);

  }


  /// Write an image view to disk.  If you supply a filename with an
  /// asterisk ('*'), each plane of the image will be saved as a
  /// separate file with the asterisk replaced with the plane number.
  template <class ImageT>
  void write_image( const std::string &filename, ImageViewBase<ImageT> const& out_image,
                    const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    ImageFormat out_image_format = out_image.format();

    unsigned files = 1;
    // If there's an asterisk, save one file per plane
    if( boost::find_last(filename,"*") ) {
      files = out_image_format.planes;
      out_image_format.planes = 1;
    }

    for( unsigned p=0; p<files; ++p ) {
      std::string name = filename;
      if( files > 1 ) boost::replace_last( name, "*",  str( boost::format("%1%") % p ) );
      VW_OUT(InfoMessage, "fileio") << "\tSaving image: " << name << "\t";
      DiskImageResource *r = DiskImageResource::create( name, out_image_format );
      VW_OUT(InfoMessage, "fileio") << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

      if ( files > 1 ) {
        write_image(*r, select_plane(out_image.impl(),p), SubProgressCallback(progress_callback,float(p)/float(files),float(p+1)/float(files)));
        progress_callback.report_finished();
      } else {
        write_image(*r, select_plane(out_image.impl(),p), progress_callback);
      }

      delete r;
    }
  }

  /// Write a std::vector of image views.  Supply a filename with an
  /// asterisk ('*') and, each image in the vector will be saved as
  /// a separate file on disk with the asterisk will be replaced with
  /// the image number.
  template <class ElemT>
  void write_image( const std::string &filename, std::vector<ElemT> const& out_image_vector ) {

    // If there's an asterisk, save one file per plane
    if( ! boost::find_last(filename,"*") ) {
      vw_throw( vw::ArgumentErr() << "write_image: filename must contain * when writing a vector of image views\n" );
    }

    for (unsigned i=0; i<out_image_vector.size(); i++){
      std::string name = filename;
      boost::replace_last( name, "*",  str( boost::format("%1%") % i ) );
      write_image( name, out_image_vector[i] );
    }
  }

  // Get the no-data value if available.
  template<class T>
  bool read_nodata_val(std::string const& file, T & nodata_val){
    boost::scoped_ptr<SrcImageResource>
      rsrc(DiskImageResource::open(file));
    if ( rsrc->has_nodata_read() ){
      nodata_val = rsrc->nodata_read();
      return true;
    }
    return false;
  }

  ImageFormat image_format(const std::string& filename);

  // Convenience function to get the width and height of an image.
  vw::Vector2i file_image_size( std::string const& input );

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOURCE_H__
