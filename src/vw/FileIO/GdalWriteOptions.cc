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
}

GdalWriteOptionsDescription::GdalWriteOptionsDescription(GdalWriteOptions& opt) {
  namespace po = boost::program_options;
  (*this).add_options()
    ("threads",      po::value(&opt.num_threads)->default_value(0),
        "Select the number of threads to use for each process. If 0, use the value in ~/.vwrc.")
    ("tile-size",  po::value(&opt.raster_tile_size)->default_value
     (Vector2i(vw_settings().default_tile_size(),
               vw_settings().default_tile_size()),"256, 256"),
        "Image tile size used for multi-threaded processing.")
    ("cache-size-mb", po::value(&opt.cache_size_mb)->default_value(1024),
        "Set the system cache size, in MB, for each process.")
    ("no-bigtiff",   "Tell GDAL to not create bigtiffs.")  // gets stored in vm.count("no-bigtiff")
    ("tif-compress", po::value(&opt.tif_compress)->default_value("LZW"),
        "TIFF Compression method. [None, LZW, Deflate, Packbits]")
    ("version,v",    "Display the version of software.")
    ("help,h",       "Display this help message.");
}
  
} // end namespace vw

