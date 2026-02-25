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

/// \file FileIO/DiskImageResource.h
///
/// An abstract base class referring to an image on disk.
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCE_H__
#define __VW_FILEIO_DISKIMAGERESOURCE_H__

#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/core/noncopyable.hpp>

#include <vw/Image/ImageResource.h>

namespace vw {

  // Return a smart pointer, this is easier to manage
  class DiskImageResource;
  boost::shared_ptr<DiskImageResource> DiskImageResourcePtr(std::string const& image_file);

  /// A base class for ImageResource objects where the buffer is on disk.
  /// - Base class from which specific file handlers derive.
  /// - Noncopyable because every impl is noncopyable
  class DiskImageResource: public ImageResource,
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
    static DiskImageResource* open(std::string const& filename);

    /// Create a new DiskImageResource of the appropriate type
    /// pointing to a newly-created empty file on disk.
    static DiskImageResource* create(std::string const& filename,
                                     ImageFormat const& format);
    static DiskImageResource* create(std::string const& filename,
                                     ImageFormat const& format,
                                     std::string const& type);

    typedef DiskImageResource* (*construct_open_func)(std::string const& filename);

    typedef DiskImageResource* (*construct_create_func)(std::string const& filename,
                                                        ImageFormat const& format);

    /// Specify that a certain file extension should be opened with a certain
    /// DiskImageResource derived type.
    static void register_file_type(std::string const& extension,
                                   std::string const& disk_image_resource_type,
                                   construct_open_func   open_func,
                                   construct_create_func create_func);

    // Specify whether values should be rescaled when converting
    // from one channel type to another during reads or writes.
    // Defaults to the value of default_rescale, which defaults to true.
    void set_rescale(bool rescale);

    // Specify a global default for rescale. This is a dangerous
    // function to use, because some higher-level Vision Workbench
    // features rely on rescaling. Use at your own risk.
    static void set_default_rescale(bool rescale);

    // TODO: This has always been the default, but it probably shouldn't be.
    virtual void flush() {}

  protected:
    DiskImageResource(std::string const& filename):
      m_filename(filename), m_rescale(default_rescale) {}
    ImageFormat m_format;
    std::string m_filename;
    bool        m_rescale;
    static bool default_rescale;
  };

  ImageFormat image_format(const std::string& filename);

  // Convenience function to get the width and height of an image.
  vw::Vector2i file_image_size(std::string const& input);

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOURCE_H__
