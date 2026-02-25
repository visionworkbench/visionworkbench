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

/// \file FileIO/DiskImageIO.h
///
/// Convenience functions for reading/writing images by filename.
/// These combine DiskImageResource (file handle) with ImageIO
/// (pixel transfer) into simple one-call operations.
///
#ifndef __VW_FILEIO_DISKIMAGEIO_H__
#define __VW_FILEIO_DISKIMAGEIO_H__

#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/scoped_ptr.hpp>

#include <vw/Core/Log.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageChannels.h>

namespace vw {

  /// Read an image on disk into a vw::ImageView<T> object.
  template <class PixelT>
  void read_image(ImageView<PixelT>& in_image, const std::string &filename) {

    VW_OUT(InfoMessage, "fileio") << "\tLoading image: " << filename << "\t";

    // Open the file for reading
    boost::shared_ptr<DiskImageResource> r(DiskImageResourcePtr(filename));

    VW_OUT(InfoMessage, "fileio") << r->cols() << "x" << r->rows()
      << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

    // Read the data
    read_image(in_image, *r);
  }

  /// Write an image view to disk. If you supply a filename with an
  /// asterisk ('*'), each plane of the image will be saved as a
  /// separate file with the asterisk replaced with the plane number.
  template <class ImageT>
  void write_image(const std::string &filename,
                   ImageViewBase<ImageT> const& out_image,
                   const ProgressCallback &progress_callback =
                     ProgressCallback::dummy_instance()) {

    ImageFormat out_image_format = out_image.format();

    unsigned files = 1;
    // If there's an asterisk, save one file per plane
    if (boost::find_last(filename, "*")) {
      files = out_image_format.planes;
      out_image_format.planes = 1;
    }

    for (unsigned p = 0; p < files; p++) {
      std::string name = filename;
      if (files > 1)
        boost::replace_last(name, "*",
                            str(boost::format("%1%") % p));
      VW_OUT(InfoMessage, "fileio") << "\tSaving image: " << name << "\t";
      DiskImageResource *r =
        DiskImageResource::create(name, out_image_format);
      VW_OUT(InfoMessage, "fileio") << r->cols() << "x" << r->rows()
        << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

      if (files > 1) {
        write_image(*r, select_plane(out_image.impl(), p),
                    SubProgressCallback(progress_callback,
                                        float(p) / float(files),
                                        float(p + 1) / float(files)));
        progress_callback.report_finished();
      } else {
        write_image(*r, select_plane(out_image.impl(), p),
                    progress_callback);
      }

      delete r;
    }
  }

  /// Write a std::vector of image views. Supply a filename with an
  /// asterisk ('*') and each image in the vector will be saved as
  /// a separate file on disk with the asterisk replaced with the
  /// image number.
  template <class ElemT>
  void write_image(const std::string &filename,
                   std::vector<ElemT> const& out_image_vector) {

    // If there's an asterisk, save one file per plane
    if (!boost::find_last(filename, "*"))
      vw_throw(vw::ArgumentErr()
        << "write_image: filename must contain * "
        << "when writing a vector of image views\n");

    for (unsigned i = 0; i < out_image_vector.size(); i++) {
      std::string name = filename;
      boost::replace_last(name, "*",
                          str(boost::format("%1%") % i));
      write_image(name, out_image_vector[i]);
    }
  }

  // Get the no-data value if available.
  template<class T>
  bool read_nodata_val(std::string const& file, T & nodata_val) {
    boost::scoped_ptr<SrcImageResource>
      rsrc(DiskImageResource::open(file));
    if (rsrc->has_nodata_read()) {
      nodata_val = rsrc->nodata_read();
      return true;
    }
    return false;
  }

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGEIO_H__
