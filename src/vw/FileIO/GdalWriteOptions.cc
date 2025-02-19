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
#include <boost/algorithm/string.hpp>

namespace boost { namespace program_options {
  // Custom value semantics, these explain how many tokens should be ingested.

  template <class T, class charT = char>
  class typed_2_value : public typed_value<T,charT> {
  public:
    typed_2_value(T* store_to) : typed_value<T,charT>(store_to) {
      this->multitoken();
    }
    
    unsigned min_tokens() const { return 2; }
    unsigned max_tokens() const { return 2; }
  };
  
  typed_2_value<vw::Vector2i>*
  value( vw::Vector2i* v ) {
    typed_2_value<vw::Vector2i>* r =
      new typed_2_value<vw::Vector2i>(v);
    return r;
  }

  // Validator for Vector2i
  template <>
  void validate(boost::any& v,
                const std::vector<std::string>& values,
                vw::Vector2i*, long ) {
    validators::check_first_occurrence(v);

    // Concatenate and then split again, so that the user can mix
    // comma and space delimited values.
    std::string joined = boost::algorithm::join(values, " ");
    std::vector<std::string> cvalues;
    boost::split(cvalues, joined, is_any_of(", "), boost::token_compress_on);

    if ( cvalues.size() != 2 )
      boost::throw_exception(invalid_syntax(invalid_syntax::missing_parameter));

    try {
      vw::Vector2i output(boost::lexical_cast<vw::int32>( cvalues[0]),
                          boost::lexical_cast<vw::int32>( cvalues[1]));
      v = output;
    } catch (boost::bad_lexical_cast const& e ) {
      boost::throw_exception(validation_error(validation_error::invalid_option_value));
    }
  }
  
}}

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

GdalWriteOptionsDescription::GdalWriteOptionsDescription(GdalWriteOptions& opt,
                                                         bool adjust_tile_size_opt) {
  namespace po = boost::program_options;
  (*this).add_options()
    ("threads",      po::value(&opt.num_threads)->default_value(0),
     "Select the number of threads to use for each process. If 0, use the value in ~/.vwrc.");
  
  if (!adjust_tile_size_opt) {
    (*this).add_options()
      ("tile-size",  po::value(&opt.raster_tile_size)->default_value
       (vw::Vector2i(vw_settings().default_tile_size(),
                     vw_settings().default_tile_size()), "256 256"), 
       "Image tile size used for multi-threaded processing.");
  } else {
    // This is needed for dem_mosaic, which already has a --tile-size option
    (*this).add_options()
      ("tif-tile-size",  po::value(&opt.raster_tile_size)->default_value
       (vw::Vector2i(vw_settings().default_tile_size(),
                     vw_settings().default_tile_size()), "256 256"), 
       "The dimensions of each block in the output image.");
  }
  (*this).add_options()
    ("cache-size-mb", po::value(&opt.cache_size_mb)->default_value(1024),
     "Set the system cache size, in MB, for each process.")
    ("no-bigtiff",   "Tell GDAL to not create bigtiffs.")  // gets stored in vm.count("no-bigtiff")
    ("tif-compress", po::value(&opt.tif_compress)->default_value("LZW"),
     "TIFF Compression method. [None, LZW, Deflate, Packbits]")
    ("version,v",    "Display the version of software.")
    ("help,h",       "Display this help message.");
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

