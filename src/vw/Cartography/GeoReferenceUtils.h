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

#ifndef __VW_CARTOGRAPHY_GEOREFERENCE_UTILS_H__
#define __VW_CARTOGRAPHY_GEOREFERENCE_UTILS_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/Cartography/Datum.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/Exception.h>

#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/Image/Manipulation.h>

/// \file GeoReferenceUtils.h Tools for working with GeoReference objects

namespace vw {
namespace cartography {

  /// Generates a new georeference which covers a sub-region of this georeference object.
  /// - Input coordinates are pixels in the corresponding image
  GeoReference crop( GeoReference const& input,
                     double upper_left_x, double upper_left_y,
                     double width=0, double height=0 );

  /// Overload of crop() that takes a bounding box object (still in pixels)
  GeoReference crop( GeoReference const& input, BBox2 const& bbox );

  /// Modify the scale in the projected coordinate to pixel coordinate transform.
  /// - A larger scale increases the number of pixels.
  GeoReference resample( GeoReference const& input, double scale_x, double scale_y );
  GeoReference resample( GeoReference const& input, double scale );

  /// Returns the distance along the Earth's surface in meters between two points.
  /// - The distance is returned in meters.
  /// - Accuracy is not perfect but should be pretty good.
  /// - Inputs should be in lon/lat degrees.
  double haversine_circle_distance(Vector2 a, Vector2 b, double radius);

  /// Estimates meters per pixel for an image.
  double get_image_meters_per_pixel(int width, int height, GeoReference const& georef);

  //---------------------------------------------------------------------------
  // Functions for writing GDAL images - multithreaded.

  template <class ImageT>
  DiskImageResourceGDAL*
  build_gdal_rsrc( const std::string &filename,
                   ImageViewBase<ImageT> const& image,
                   GdalWriteOptions const& opt) {
    return new DiskImageResourceGDAL(filename, image.impl().format(),
                                         opt.raster_tile_size, opt.gdal_options);
  }

  /// Multi-threaded block write image with, if available, nodata, georef, and
  /// keywords to geoheader.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              bool has_georef,
                              cartography::GeoReference const& georef,
                              bool has_nodata, double nodata,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback =
                              ProgressCallback::dummy_instance(),
                              std::map<std::string, std::string> const& keywords =
                              std::map<std::string, std::string>() );

  /// Block write image without georef and nodata.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback =
                              ProgressCallback::dummy_instance(),
                              std::map<std::string, std::string> const& keywords =
                              std::map<std::string, std::string>() );

  /// Block write image with nodata.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              double nodata,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback =
                              ProgressCallback::dummy_instance(),
                              std::map<std::string, std::string> const& keywords =
                              std::map<std::string, std::string>() );

  //---------------------------------------------------------------------------
  // Functions for writing GDAL images - single threaded.

  /// Single-threaded write image with, if available, nodata, georef, and
  /// keywords to geoheader.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        bool has_georef,
                        cartography::GeoReference const& georef,
                        bool has_nodata, double nodata,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback
                        = ProgressCallback::dummy_instance(),
                        std::map<std::string, std::string> const& keywords
                        = std::map<std::string, std::string>() );


  /// Single-threaded write image with georef and keywords to geoheader.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        cartography::GeoReference const& georef,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback
                        = ProgressCallback::dummy_instance(),
                        std::map<std::string, std::string> const& keywords
                        = std::map<std::string, std::string>() );

  /// Single-threaded write image with georef, nodata, and keywords to geoheader.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        cartography::GeoReference const& georef,
                        double nodata,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback
                        = ProgressCallback::dummy_instance(),
                        std::map<std::string, std::string> const& keywords
                        = std::map<std::string, std::string>() );

  /// Single-threaded write image.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback
                        = ProgressCallback::dummy_instance(),
                        std::map<std::string, std::string> const& keywords
                        = std::map<std::string, std::string>() );


//=============================================================================
//=============================================================================
// Function definitions

  // Multi-threaded block write image with, if available, nodata, georef, and
  // keywords to geoheader.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              bool has_georef,
                              cartography::GeoReference const& georef,
                              bool has_nodata, double nodata,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback,
                              std::map<std::string, std::string> const& keywords) {

    boost::scoped_ptr<DiskImageResourceGDAL>
      rsrc(build_gdal_rsrc(filename, image, opt));

    if (has_nodata)
      rsrc->set_nodata_write(nodata);

    for (std::map<std::string, std::string>::const_iterator i = keywords.begin();
         i != keywords.end(); i++)
      cartography::write_header_string(*rsrc, i->first, i->second);

    if (has_georef)
      cartography::write_georeference(*rsrc, georef);

    block_write_image(*rsrc, image.impl(), progress_callback, opt.num_threads);
  }

  // Block write image without georef and nodata.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback,
                              std::map<std::string, std::string> const& keywords) {
    bool has_nodata = false;
    bool has_georef = false;
    float nodata = std::numeric_limits<float>::quiet_NaN();
    cartography::GeoReference georef;
    block_write_gdal_image(filename, image, has_georef, georef,
                           has_nodata, nodata, opt,
                           progress_callback, keywords);
  }


  // Block write image with nodata.
  template <class ImageT>
  void block_write_gdal_image(const std::string &filename,
                              ImageViewBase<ImageT> const& image,
                              double nodata,
                              GdalWriteOptions const& opt,
                              ProgressCallback const& progress_callback,
                              std::map<std::string, std::string> const& keywords) {
    bool has_nodata = true;
    bool has_georef = false;
    cartography::GeoReference georef;
    block_write_gdal_image(filename, image, has_georef, georef,
                           has_nodata, nodata, opt,
                           progress_callback, keywords);
  }

  // Single-threaded write functions.

  // Single-threaded write image with, if available, nodata, georef, and
  // keywords to geoheader.

  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        bool has_georef,
                        cartography::GeoReference const& georef,
                        bool has_nodata, double nodata,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback,
                        std::map<std::string, std::string> const& keywords) {

    boost::scoped_ptr<DiskImageResourceGDAL>
      rsrc( build_gdal_rsrc( filename, image, opt ) );

    if (has_nodata)
      rsrc->set_nodata_write(nodata);

    for (std::map<std::string, std::string>::const_iterator i = keywords.begin();
         i != keywords.end(); i++)
      cartography::write_header_string(*rsrc, i->first, i->second);

    if (has_georef)
      cartography::write_georeference(*rsrc, georef);

    write_image( *rsrc, image.impl(), progress_callback );
  }


  // Single-threaded write image with georef and keywords to geoheader.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        cartography::GeoReference const& georef,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback,
                        std::map<std::string, std::string> const& keywords) {

    bool has_georef = true;
    bool has_nodata = false;
    float nodata = std::numeric_limits<float>::quiet_NaN();

    write_gdal_image(filename, image, has_georef, georef,
                     has_nodata, nodata, opt, progress_callback,
                     keywords);
  }

  // Single-threaded write image with georef, nodata, and keywords to geoheader.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        cartography::GeoReference const& georef,
                        double nodata,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback,
                        std::map<std::string, std::string> const& keywords) {

    bool has_georef = true;
    bool has_nodata = true;
    write_gdal_image(filename, image, has_georef, georef,
                     has_nodata, nodata, opt, progress_callback,
                     keywords);
  }

  // Single-threaded write image.
  template <class ImageT>
  void write_gdal_image(const std::string &filename,
                        ImageViewBase<ImageT> const& image,
                        GdalWriteOptions const& opt,
                        ProgressCallback const& progress_callback,
                        std::map<std::string, std::string> const& keywords) {

    bool has_nodata = false;
    bool has_georef = false;
    float nodata = std::numeric_limits<float>::quiet_NaN();
    cartography::GeoReference georef;
    write_gdal_image(filename, image, has_georef, georef,
                     has_nodata, nodata, opt, progress_callback,
                     keywords);
  }

  /// The following namespace contains functions that return GeoReferences
  /// for certain well-known output styles, such as KML (and related
  /// functions involved in doing so).
  namespace output {
    namespace kml {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        Vector2 pos      = tx.forward( pixel );
        Vector2 x_vector = tx.forward( pixel+Vector2(1,0) ) - pos;
        Vector2 y_vector = tx.forward( pixel+Vector2(0,1) ) - pos;
        double degrees_per_pixel        = (std::min)( norm_2(x_vector), norm_2(y_vector) );
        double pixels_per_circumference = 360.0 / degrees_per_pixel;
        int scale_exponent = (int) ceil( log(pixels_per_circumference)/log(2.0) );
        if (scale_exponent >= 31) scale_exponent = 30;
        return 1 << scale_exponent;
      }
    } // namespace: vw::cartography::output::kml

    namespace tms {
      // Returns the number of pixels per planetary circumference,
      // rounding up to a power of two.
      template <class TransformT>
      inline int32 compute_resolution( TransformT const& tx, Vector2 const& pixel ) {
        // It's exactly the same as the one for KML.
        return vw::cartography::output::kml::compute_resolution(tx, pixel);
      }
    } // namespace vw::cartography::output::tms
  } // namespace vw::cartography::output

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCE_UTILS_H__
