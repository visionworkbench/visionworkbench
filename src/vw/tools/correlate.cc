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

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Debugging.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/Math/Vector.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Transform.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/MatcherIO.h>
#include <vw/Stereo/CorrelationView.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::stereo;

// TODO(oalexan1): Make all options below use the Options structure
struct Options: vw::GdalWriteOptions {};

int main( int argc, char *argv[] ) {

  Options opt;
  std::string left_file_name, right_file_name;
  float log;
  int32 h_corr_min, h_corr_max;
  int32 v_corr_min, v_corr_max;
  int32 xkernel, ykernel;
  int   lrthresh;
  int   min_lr_level;
  int   nThreads;
  int   tile_size;
  int   collar_size;
  int   sgm_filter_size;
  int   max_pyramid_levels;
  int   cost_mode;
  int   stereo_algorithm;
  bool  found_alignment = false;
  int   blob_filter_area;
  Matrix3x3 alignment;
  float mask_value;
  int   filter_radius;
  int   mem_limit_gb;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("left",       po::value(&left_file_name),                 "Explicitly specify the \"left\" input file")
    ("right",      po::value(&right_file_name),                "Explicitly specify the \"right\" input file")
    ("log",        po::value(&log)->default_value(1.4),        "Apply LOG filter with the given sigma, or 0 to disable")
    ("h-corr-min", po::value(&h_corr_min)->default_value(-30), "Minimum horizontal disparity")
    ("h-corr-max", po::value(&h_corr_max)->default_value( 30), "Maximum horizontal disparity")
    ("v-corr-min", po::value(&v_corr_min)->default_value(-5),  "Minimum vertical disparity")
    ("v-corr-max", po::value(&v_corr_max)->default_value(5),   "Maximum vertical disparity")
    ("xkernel",    po::value(&xkernel)->default_value(15),     "Horizontal correlation kernel size")
    ("ykernel",    po::value(&ykernel)->default_value(15),     "Vertical correlation kernel size")
    ("lrthresh",   po::value(&lrthresh)->default_value(2),     "Left/right correspondence threshold")
    ("min-lr-level",    po::value(&min_lr_level)->default_value(1),     "Min level to check L/R correspondence at (SGM only).")
    ("filter-radius",      po::value(&filter_radius)->default_value(5), "Disp filtering radius")
    ("cost-mode",          po::value(&cost_mode)->default_value(0), 
                        "0 - Abs difference; 1 - Sq Difference; 2 - NormXCorr; 3 - Census; 4 - Ternary Census")
    ("stereo-algorithm",   po::value(&stereo_algorithm)->default_value(0), 
                    "Choose the stereo algorithm: 0=BM, 1=SGM, 2=MGM, 3=FinalMGM")
    //("affine-subpix", "Enable affine adaptive sub-pixel correlation (slower, but more accurate)") // TODO: Unused!
    ("blob-filter-area",   po::value(&blob_filter_area)->default_value(0),     "Filter blobs of this size.")
    ("threads",            po::value(&nThreads)->default_value(0),    "Manually specify the number of threads")
    ("tile-size",          po::value(&tile_size)->default_value(0),   "Manually specify the tile size")
    ("collar-size",        po::value(&collar_size)->default_value(0), "Specify a collar size")
    ("max-mem-GB",         po::value(&mem_limit_gb)->default_value(0), "Specify the maximum memory size")
    ("sgm-filter-size",    po::value(&sgm_filter_size)->default_value(0), "Filter SGM subpixel results with this size")
    ("mask-value",         po::value(&mask_value)->default_value(-32768), "Specify a mask value")
    ("max-pyramid-levels", po::value(&max_pyramid_levels)->default_value(5),
      "Limit the maximum number of pyramid levels")
    ("debug",      "Write out debugging images")
    ;

  general_options.add(vw::GdalWriteOptionsDescription(opt));
  
  po::positional_options_description p;
  p.add("left", 1);
  p.add("right", 1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(general_options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occurred while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << general_options << std::endl;
    return 1;
  }

  if (vm.count("help")) {
    vw_out() << general_options << std::endl;
    return 1;
  }

  opt.setVwSettingsFromOpt();

  if( vm.count("left") != 1 || vm.count("right") != 1 ) {
    vw_out() << "Error: Must specify one (and only one) left and right input file!" << std::endl;
    vw_out() << general_options << std::endl;
    return 1;
  }

  if (nThreads > 0)
    vw::vw_settings().set_default_num_threads(nThreads);
  if (tile_size > 0)
    vw::vw_settings().set_default_tile_size(tile_size);

  // If the user provided a match file, use it to align the images.
  std::string match_file = vw::ip::match_filename("", left_file_name, right_file_name);
  if (fs::exists(match_file)) {
    vw_out() << "Found a match file. Using it to pre-align images.\n";
    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
    ip::read_binary_match_file(match_file, matched_ip1, matched_ip2);
    std::vector<Vector3> ransac_ip1 = ip::iplist_to_vectorlist(matched_ip1);
    std::vector<Vector3> ransac_ip2 = ip::iplist_to_vectorlist(matched_ip2);
    vw::math::RandomSampleConsensus<vw::math::HomographyFittingFunctor, vw::math::InterestPointErrorMetric> 
        ransac( vw::math::HomographyFittingFunctor(), vw::math::InterestPointErrorMetric(), 100, 30, ransac_ip1.size()/2, true );
    alignment = ransac( ransac_ip2, ransac_ip1 );

    DiskImageView<PixelGray<float> > right_disk_image( right_file_name );
    right_file_name = "aligned_right.tif";
    write_image( right_file_name, transform(right_disk_image, HomographyTransform(alignment)),
                 TerminalProgressCallback( "tools.correlate", "Aligning: ") );
    found_alignment = true;
  }

  // Load input images
  DiskImageView<PixelGray<float> > left_disk_image (left_file_name );
  DiskImageView<PixelGray<float> > right_disk_image(right_file_name );
  int cols = std::min(left_disk_image.cols(),right_disk_image.cols());
  int rows = std::min(left_disk_image.rows(),right_disk_image.rows());
  ImageViewRef<PixelGray<float> > left  = edge_extend(left_disk_image, 0,0,cols,rows);
  ImageViewRef<PixelGray<float> > right = edge_extend(right_disk_image,0,0,cols,rows);

  // Set up masks
  ImageView<uint8> left_mask  = constant_view( uint8(255), left );
  ImageView<uint8> right_mask = constant_view( uint8(255), right);
  if (vm.count("mask-value")) {
    left_mask  = apply_mask(copy_mask(left_mask,  create_mask(left,  mask_value)), 0);
    right_mask = apply_mask(copy_mask(right_mask, create_mask(right, mask_value)), 0);
    
    write_image("correlate_left_mask.tif", left_mask);
    write_image("correlate_right_mask.tif", right_mask);
  }

  bool write_debug_images = (vm.count("debug"));

  // TODO: Hook up to options!
  SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = SemiGlobalMatcher::SUBPIXEL_LC_BLEND;
  Vector2i sgm_search_buffer(4,2);
  size_t memory_limit_mb = 1024*5;
  if (mem_limit_gb > 0)
    memory_limit_mb = mem_limit_gb * 1024;

  int corr_timeout = 0;
  double seconds_per_op = 0.0;
  BBox2i search_range(Vector2i(h_corr_min, v_corr_min), 
                      Vector2i(h_corr_max, v_corr_max));
  Vector2i kernel_size(xkernel, ykernel);
  std::cout << "Correlate max search range = " << search_range << std::endl;

  // Here we don't record the L-R disparity difference
  ImageView<PixelMask<float>> * lr_disp_diff = NULL;
  Vector2i region_ul = Vector2i(0, 0);
  
  ImageViewRef<PixelMask<Vector2f>> disparity_map
    = stereo::pyramid_correlate(left, right,
                                left_mask, right_mask,
                                stereo::PREFILTER_LOG, log,
                                search_range, kernel_size,
                                static_cast<vw::stereo::CostFunctionType>(cost_mode), 
                                corr_timeout, seconds_per_op,
                                lrthresh, min_lr_level,
                                filter_radius, max_pyramid_levels, 
                                static_cast<vw::stereo::CorrelationAlgorithm>(stereo_algorithm), 
                                collar_size,
                                sgm_subpixel_mode,
                                sgm_search_buffer,
                                memory_limit_mb,
                                blob_filter_area,
                                lr_disp_diff, region_ul, 
                                write_debug_images);
  
  // TODO: Call the code used in stereo_fltr!

/*

// TODO: Debug code for experimenting with disparity filtering code.
  int texture_smooth_range = 15;
  float texture_max_percentage = 0.85;
  int max_kernel_size = 13;

  std::cout << "Generating texture image...\n";
  ImageView<float> texture_image;
  ImageView<float> leftR = left;
  float max_texture = texture_measure(leftR, texture_image, texture_smooth_range);
  write_image( "texture_image.tif", texture_image );

  std::cout << "Rasterizing disparity image...\n";
  ImageView<PixelMask<Vector2f> > disparity_map_raster = disparity_map;
  //disparity_map.rasterize(disparity_map_raster, bounding_box(disparity_map));
  
  std::cout << "Filtering disparity image...\n";
  ImageView<PixelMask<Vector2f> > disparity_map_filtered;
  texture_preserving_disparity_filter(disparity_map_raster, disparity_map_filtered, texture_image, 
                                      texture_max_percentage*max_texture, max_kernel_size);
  
  write_image( "texture_filtered_disp.tif", disparity_map_filtered );
  //write_image( "subpixel_disp.tif", disparity );
  //write_image( "subpixel_deltas.tif", deltas );
  
  std::cout << "Done!\n";


*/

  //ImageViewRef<PixelMask<Vector2f> > result = pixel_cast<PixelMask<Vector2f> >(disparity_map);
  if ( found_alignment )
    disparity_map = transform_disparities(disparity_map, HomographyTransform(alignment) );

  // Actually invoke the raster
  {
    vw::Timer corr_timer("Correlation Time");
    vw::GdalWriteOptions geo_opt;
    geo_opt.raster_tile_size = Vector2i(1024, 1024);
    if (stereo_algorithm == VW_CORRELATION_BM) {
      vw::cartography::block_write_gdal_image("disparity.tif", disparity_map, geo_opt);
    }
    else { // SGM/MGM needs to be rasterized in a single tile.
      ImageView<PixelMask<Vector2f> > result = disparity_map;
      vw::cartography::block_write_gdal_image("disparity.tif", result, geo_opt);
    }
  }

  //// Write disparity debug images
  //DiskImageView<PixelMask<Vector2i> > solution("disparity.tif");
  //BBox2 disp_range = get_disparity_range(solution);
  //std::cout << "Found disparity range: " << disp_range << "\n";

  // Are these working properly?
  //write_image( "x_disparity.tif",
  //             channel_cast<uint8>(apply_mask(copy_mask(clamp(normalize(select_channel(solution,0), 
  //                                 disp_range.min().x(), disp_range.max().x(),0,255)),solution))) );
  //write_image( "y_disparity.tif",
  //             channel_cast<uint8>(apply_mask(copy_mask(clamp(normalize(select_channel(solution,1), 
  //                                 disp_range.min().y(), disp_range.max().y(),0,255)),solution))) );

  return 0;
}
