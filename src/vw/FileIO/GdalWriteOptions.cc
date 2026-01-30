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

#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <boost/algorithm/string.hpp>

namespace vw {

GdalWriteOptions::GdalWriteOptions() {
#if defined(VW_HAS_BIGTIFF) && VW_HAS_BIGTIFF == 1
  gdal_options["BIGTIFF"]  = "YES";
#else
  gdal_options["BIGTIFF"] = "NO";
#endif
  gdal_options["COMPRESS"] = "LZW";

  // The values below will be overwritten by the boost options
  raster_tile_size = Vector2i(vw_settings().default_tile_size(),
                              vw_settings().default_tile_size());
  num_threads      = vw_settings().default_num_threads();
  cache_size_mb    = vw_settings().system_cache_size() / (1024.0 * 1024.0); // bytes to MB
  cog              = false;
}

void GdalWriteOptions::setVwSettingsFromOpt() {
    
  // If the user did not set the number of threads, use what is set in
  // .vwrc.
  if (this->num_threads <= 0)
    this->num_threads = vw_settings().default_num_threads();
  
  // Same for cache size
  double MB = 1024.0 * 1024.0; // 1 MB in bytes
  if (this->cache_size_mb <= 0)
    this->cache_size_mb = vw_settings().system_cache_size() / MB;
  
  // Here we ensure that this->num_threads and default_num_threads()
  // are consistent among themselves. Same for cache size.
  vw::vw_settings().set_default_num_threads(this->num_threads);
  vw::vw_settings().set_system_cache_size(this->cache_size_mb * MB);
  
  // Print the message below just once per process.
  static bool verbose = true;
  if (verbose) {
    vw::vw_out() << "\t--> Setting number of processing threads to: "
                 <<  vw_settings().default_num_threads() << std::endl;
    verbose = false;
  }
  
  boost::algorithm::to_upper(this->tif_compress);
  boost::algorithm::trim( this->tif_compress );
  VW_ASSERT( this->tif_compress == "NONE" || this->tif_compress == "LZW" ||
             this->tif_compress == "DEFLATE" || this->tif_compress == "PACKBITS",
             ArgumentErr() << "\"" << this->tif_compress
             << "\" is not a valid options for TIF_COMPRESS." );
  this->gdal_options["COMPRESS"] = this->tif_compress;
  
}
  
} // end namespace vw

