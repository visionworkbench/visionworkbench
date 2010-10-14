// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <vw/Core/Debugging.h>
#include <vw/Math.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Stereo/ReferenceCorrelator.h>
#include <vw/Stereo/PyramidCorrelator.h>

using namespace vw;
using namespace vw::stereo;

int main( int argc, char *argv[] ) {
  try {

    std::string left_file_name, right_file_name;
    float log, slog;
    int h_corr_min, h_corr_max;
    int v_corr_min, v_corr_max;
    int xkernel, ykernel;
    int lrthresh;
    float corrscore_thresh;
    int cost_blur;
    int correlator_type;
    bool found_alignment = false;
    Matrix3x3 alignment;

    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Display this help message")
      ("left", po::value(&left_file_name), "Explicitly specify the \"left\" input file")
      ("right", po::value(&right_file_name), "Explicitly specify the \"right\" input file")
      ("slog", po::value(&slog)->default_value(1.0), "Apply SLOG filter with the given sigma, or 0 to disable")
      ("log", po::value(&log)->default_value(0.0), "Apply LOG filter with the given sigma, or 0 to disable")
      ("h-corr-min", po::value(&h_corr_min)->default_value(0), "Minimum horizontal disparity")
      ("h-corr-max", po::value(&h_corr_max)->default_value(0), "Maximum horizontal disparity")
      ("v-corr-min", po::value(&v_corr_min)->default_value(5), "Minimum vertical disparity")
      ("v-corr-max", po::value(&v_corr_max)->default_value(5), "Maximum vertical disparity")
      ("xkernel", po::value(&xkernel)->default_value(15), "Horizontal correlation kernel size")
      ("ykernel", po::value(&ykernel)->default_value(15), "Vertical correlation kernel size")
      ("lrthresh", po::value(&lrthresh)->default_value(2), "Left/right correspondence threshold")
      ("csthresh", po::value(&corrscore_thresh)->default_value(1.0), "Correlation score rejection threshold (1.0 is Off <--> 2.0 is Aggressive outlier rejection")
      ("cost-blur", po::value(&cost_blur)->default_value(1), "Kernel size for bluring the cost image")
      ("correlator-type", po::value(&correlator_type)->default_value(0), "0 - Abs difference; 1 - Sq Difference; 2 - NormXCorr")
      ("hsubpix", "Enable horizontal sub-pixel correlation")
      ("vsubpix", "Enable vertical sub-pixel correlation")
      ("affine-subpix", "Enable affine adaptive sub-pixel correlation (slower, but more accurate)")
      ("reference", "Use the slower, simpler reference correlator")
      ("pyramid", "Use the pyramid based correlator")
      ("bitimage", "Force the use of the optimized bit-image correlator")
      ("nonbitimage", "Fore the use of the slower, non bit-image optimized correlator")
      ;
    po::positional_options_description p;
    p.add("left", 1);
    p.add("right", 1);

    po::variables_map vm;
    try {
      po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
      po::notify( vm );
    } catch (po::error &e) {
      std::cout << "An error occured while parsing command line arguments.\n";
      std::cout << "\t" << e.what() << "\n\n";
      std::cout << desc << std::endl;
      return 1;
    }

    if( vm.count("help") ) {
      vw_out() << desc << std::endl;
      return 1;
    }

    if( vm.count("left") != 1 || vm.count("right") != 1 ) {
      vw_out() << "Error: Must specify one (and only one) left and right input file!" << std::endl;
      vw_out() << desc << std::endl;
      return 1;
    }

    std::string match_filename = fs::path( left_file_name ).replace_extension().string() + "__" + fs::path( right_file_name ).stem() + ".match";
    if ( fs::exists( match_filename ) ) {
      vw_out() << "Found a match file. Using it to pre-align images.\n";
      std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
      ip::read_binary_match_file( match_filename,
                                  matched_ip1, matched_ip2 );
      std::vector<Vector3> ransac_ip1 = ip::iplist_to_vectorlist(matched_ip1);
      std::vector<Vector3> ransac_ip2 = ip::iplist_to_vectorlist(matched_ip2);
      vw::math::RandomSampleConsensus<vw::math::HomographyFittingFunctor, vw::math::InterestPointErrorMetric> ransac( vw::math::HomographyFittingFunctor(), vw::math::InterestPointErrorMetric(), 30 );
      alignment = ransac( ransac_ip2, ransac_ip1 );

      DiskImageView<PixelGray<float> > right_disk_image( right_file_name );
      right_file_name = "aligned_right.tif";
      write_image( right_file_name, transform(right_disk_image, HomographyTransform(alignment)),
                   TerminalProgressCallback( "tools.correlate", "Aligning: ") );
      found_alignment = true;
    }

    DiskImageView<PixelGray<float> > left_disk_image(left_file_name );
    DiskImageView<PixelGray<float> > right_disk_image(right_file_name );
    int cols = std::min(left_disk_image.cols(),right_disk_image.cols());
    int rows = std::min(left_disk_image.rows(),right_disk_image.rows());
    ImageViewRef<PixelGray<float> > left = edge_extend(left_disk_image,0,0,cols,rows);
    ImageViewRef<PixelGray<float> > right = edge_extend(right_disk_image,0,0,cols,rows);

    ImageViewRef<uint8> left_mask = select_channel(edge_mask(pixel_cast_rescale<uint8>(left)),1);
    ImageViewRef<uint8> right_mask = select_channel(edge_mask(pixel_cast_rescale<uint8>(right)),1);

    stereo::CorrelatorType corr_type = ABS_DIFF_CORRELATOR; // the default
    if (correlator_type == 1)
      corr_type = SQR_DIFF_CORRELATOR;
    else if (correlator_type == 2)
      corr_type = NORM_XCORR_CORRELATOR;

    ImageView<PixelMask<Vector2f> > disparity_map;
    if (vm.count("reference")) {
      vw::stereo::ReferenceCorrelator correlator( h_corr_min, h_corr_max,
                                                  v_corr_min, v_corr_max,
                                                  xkernel, ykernel,
                                                  true, lrthresh,
                                                  (vm.count("hsubpix")>0),
                                                  (vm.count("vsubpix")>0),
                                                  (vm.count("affine-subpix")>0) );
      if (log > 0)
        disparity_map = correlator( left, right, stereo::LogStereoPreprocessingFilter(log));
      else
        disparity_map = correlator( left, right, stereo::SlogStereoPreprocessingFilter(slog));
    } else if (vm.count("pyramid")) {
      vw::stereo::PyramidCorrelator correlator( BBox2(Vector2(h_corr_min, v_corr_min),
                                                      Vector2(h_corr_max, v_corr_max)),
                                                Vector2i(xkernel, ykernel),
                                                lrthresh,
                                                corrscore_thresh,
                                                cost_blur,
                                                corr_type);
      correlator.set_debug_mode("debug");
      {
        vw::Timer corr_timer("Correlation Time");
        if (log > 0)
          disparity_map = correlator( left, right, left_mask, right_mask, stereo::LogStereoPreprocessingFilter(log));
        else
          disparity_map = correlator( left, right, left_mask, right_mask, stereo::SlogStereoPreprocessingFilter(slog));
      }
    } else {
      vw::stereo::OptimizedCorrelator correlator( BBox2i(Vector2(h_corr_min, v_corr_min),
                                                         Vector2(h_corr_max, v_corr_max)),
                                                  xkernel,
                                                  lrthresh,
                                                  corrscore_thresh,
                                                  cost_blur,
                                                  corr_type);
      {
        vw::Timer corr_timer("Correlation Time");
        if (log > 0)
          disparity_map = correlator( left, right, stereo::LogStereoPreprocessingFilter(log));
        else
          disparity_map = correlator( left, right, stereo::SlogStereoPreprocessingFilter(slog));
      }
    }

    ImageViewRef<PixelMask<Vector2f> > result = disparity_map;
    if ( found_alignment )
      result = transform_disparities(disparity_map, HomographyTransform(alignment) );

    // Write disparity debug images
    BBox2 disp_range = get_disparity_range(result);
    std::cout << "Found disparity range: " << disp_range << "\n";
    ImageViewRef<float32> horizontal = apply_mask(copy_mask(clamp(normalize(select_channel(result,0), disp_range.min().x(), disp_range.max().x(),0,1)),result));
    ImageViewRef<float32> vertical = apply_mask(copy_mask(clamp(normalize(select_channel(result,1), disp_range.min().y(), disp_range.max().y(),0,1)),result));

    write_image( "x_disparity.png", channel_cast_rescale<uint8>(horizontal) );
    write_image( "y_disparity.png", channel_cast_rescale<uint8>(vertical) );
  }
  catch( vw::Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}

