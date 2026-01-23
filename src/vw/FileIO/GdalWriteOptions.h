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


/// \file GdalWriteOptions.h
///
/// Boost program options to add to each tool shipped with ASP. Some are
/// passed to GDAL, such as the image tile size.

/// Advanced users can pass custom options to GDAL when creating a
/// resource.  Here is an example showing how to make a tiled,
/// compressed BigTIFF (assuming libtiff 4.0 or greater):
///
///   DiskImageResourceGDAL::Options options;
///   options["COMPRESS"] = "LZW";
///   options["BIGTIFF"] = "YES";
///   DiskImageResourceGDAL resource( "filename.tif",
///                                   image.format(),
///                                   Vector2i(256,256),
///                                   options );
///   write_image( resource, image );
///
#ifndef __VW_FILEIO_GDALWRITEOPTIONS_H__
#define __VW_FILEIO_GDALWRITEOPTIONS_H__

#include <vw/vw_config.h>
#include <string>
#include <map>

#include <vw/FileIO/DiskImageResourceGDAL.h>

namespace vw {

  /// Standard options for multi-threaded GDAL (tif) image writing.
  /// - num_threads sets the number of parallel block-writing threads when calling one
  ///   of the block write functions in this file.  By default it is set to
  ///   vw_settings().default_num_threads().
  struct GdalWriteOptions {
    DiskImageResourceGDAL::Options gdal_options;
    Vector2i     raster_tile_size;
    int32        num_threads;
    std::int64_t cache_size_mb; // Set the cache size, in MB. Modifies VW's system_cache_size.
    std::string  tif_compress;

    GdalWriteOptions();

    void setVwSettingsFromOpt();
  };
  
} // namespace vw

#endif // __VW_FILEIO_GDALWRITEOPTIONS_H__
