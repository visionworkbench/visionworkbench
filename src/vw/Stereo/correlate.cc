#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Stereo/ReferenceCorrelator.h>
#include <vw/Stereo/MultiresolutionCorrelator.h>
#include <vw/Stereo/Correlator.h>
#include <vw/Stereo/BlockCorrelator.h>

using namespace vw;
using namespace vw::stereo;
using namespace vw::stereo::disparity;

int main( int argc, char *argv[] ) {
  try {

    std::string left_file_name, right_file_name;
    float log, slog;
    int xoffset, yoffset;
    int xrange, yrange;
    int xkernel, ykernel;
    int lrthresh;
    float corrscore_thresh;
    int block_size = 0;

    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("left", po::value<std::string>(&left_file_name), "Explicitly specify the \"left\" input file")
      ("right", po::value<std::string>(&right_file_name), "Explicitly specify the \"right\" input file")
      ("slog", po::value<float>(&slog)->default_value(1.0), "Apply SLOG filter with the given sigma, or 0 to disable")
      ("log", po::value<float>(&log)->default_value(0.0), "Apply LOG filter with the given sigma, or 0 to disable")
      ("xoffset", po::value<int>(&xoffset)->default_value(0), "Overall horizontal offset between images")
      ("yoffset", po::value<int>(&yoffset)->default_value(0), "Overall vertical offset between images")
      ("xrange", po::value<int>(&xrange)->default_value(5), "Allowed range of horizontal disparity")
      ("yrange", po::value<int>(&yrange)->default_value(5), "Allowed range of vertical disparity")
      ("xkernel", po::value<int>(&xkernel)->default_value(15), "Horizontal correlation kernel size")
      ("ykernel", po::value<int>(&ykernel)->default_value(15), "Vertical correlation kernel size")
      ("lrthresh", po::value<int>(&lrthresh)->default_value(2), "Left/right correspondence threshold")
      ("csthresh", po::value<float>(&corrscore_thresh)->default_value(1.0), "Correlation score rejection threshold (1.0 is Off <--> 2.0 is Aggressive outlier rejection")
      ("block", po::value<int>(&block_size), "Use the prototype block image correlator")
      ("hsubpix", "Enable horizontal sub-pixel correlation")
      ("vsubpix", "Enable vertical sub-pixel correlation")
      ("reference", "Use the slower, simpler reference correlator")
      ("pyramid", "Use the pyramid based correlator")
      ("multiresolution", "Use the prototype multiresolution correlator")
      ("bitimage", "Force the use of the optimized bit-image correlator")
      ("nonbitimage", "Fore the use of the slower, non bit-image optimized correlator")
      ;
    po::positional_options_description p;
    p.add("left", 1);
    p.add("right", 1);

    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );

    if( vm.count("help") ) {
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( vm.count("left") != 1 || vm.count("right") != 1 ) {
      std::cout << "Error: Must specify one (and only one) left and right input file!" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }

    vw::ImageView<vw::PixelGray<float> > left, right;
    read_image( left, left_file_name );
    read_image( right, right_file_name );
    int cols = std::max(left.cols(),right.cols());
    int rows = std::max(left.rows(),right.rows());
    left = edge_extend(left,0,0,cols,rows);
    right = edge_extend(right,0,0,cols,rows);

    if( log > 0 && slog > 0 ) {
      std::cout << "Error: Can only specify one of LOG and SLOG pre-filter!" << std::endl;
    }

    bool bit_image = false;
    if( log > 0 && !vm.count("pyramid") ) {
      std::cout << "Applying LOG filter..." << std::endl;
      left = laplacian_filter( gaussian_filter( left, log ) );
      right = laplacian_filter( gaussian_filter( right, log ) );
    }

    if( slog > 0.0 && !vm.count("pyramid") ) {
      bit_image = true;
      std::cout << "Applying SLOG filter..." << std::endl;
      left = threshold( laplacian_filter( gaussian_filter( left, slog ) ) );
      right = threshold( laplacian_filter( gaussian_filter( right, slog ) ) );
    }

    ImageView<PixelDisparity<float> > disparity_map;
    if (vm.count("reference")>0) {
      vw::stereo::ReferenceCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                  yoffset-yrange, yoffset+yrange,
                                                  xkernel, ykernel,
                                                  true, lrthresh,
                                                  (vm.count("hsubpix")>0),
                                                  (vm.count("vsubpix")>0) );
      disparity_map = correlator( left, right, bit_image );
    } else if (vm.count("pyramid")>0) {
      if (slog == 0) {
        std::cout << "Pyramid correlator cannot operate with an slog of 0.   Please choose a postive value.\n";
        exit(0);
      }
      vw::stereo::Correlator correlator( BBox2(Vector2(xoffset-xrange, yoffset-yrange),
                                               Vector2(xoffset+xrange, yoffset+yrange)),
                                         Vector2i(xkernel, ykernel),
                                         slog, lrthresh,
                                         corrscore_thresh,
                                         (vm.count("hsubpix")>0),
                                         (vm.count("vsubpix")>0) );
      correlator.set_debug_mode("debug");
      vw::Timer corr_timer("Correlation Time");
      disparity_map = correlator( left, right );      
    } else if (vm.count("multiresolution")>0) {
      read_image( left, left_file_name );
      read_image( right, right_file_name );
      int cols = std::max(left.cols(),right.cols());
      int rows = std::max(left.rows(),right.rows());
      left = edge_extend(left,0,0,cols,rows);
      right = edge_extend(right,0,0,cols,rows);
      vw::stereo::MultiresolutionCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                        yoffset-yrange, yoffset+yrange,
                                                        xkernel, ykernel,
                                                        true, lrthresh, 
                                                        corrscore_thresh,
                                                        slog,
                                                        (vm.count("hsubpix")>0),
                                                        (vm.count("vsubpix")>0) );
      vw::Timer corr_timer("Correlation Time");
      disparity_map = correlator( left, right );      
    } else if (block_size != 0) {
      vw::stereo::BlockCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                              yoffset-yrange, yoffset+yrange,
                                              xkernel, ykernel,
                                              true, lrthresh, corrscore_thresh, 
                                              block_size,
                                              (vm.count("hsubpix")>0),
                                              (vm.count("vsubpix")>0) );
      if (vm.count("bitimage")) {
        std::cout << "Forcing the use of the bit-image optimized correlator.\n";
        bit_image = true;
      } else if (vm.count("nonbitimage")) {
        std::cout << "Forcing the use of the non bit-image optimized correlator.\n";
        bit_image = false;
      } 

      if (bit_image) {
        ImageView<PixelGray<uint8> > left_bitimage = channel_cast<uint8>(left);
        ImageView<PixelGray<uint8> > right_bitimage = channel_cast<uint8>(right);
        vw::Timer corr_timer("Correlation Time");
        disparity_map = correlator( left_bitimage, right_bitimage, bit_image );
      } else {
        vw::Timer corr_timer("Correlation Time");
        disparity_map = correlator( left, right, bit_image );
      }

    } else {
      vw::stereo::OptimizedCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                  yoffset-yrange, yoffset+yrange,
                                                  xkernel, ykernel,
                                                  true, lrthresh, corrscore_thresh,
                                                  (vm.count("hsubpix")>0),
                                                  (vm.count("vsubpix")>0) );
      if (vm.count("bitimage")) {
        std::cout << "Forcing the use of the bit-image optimized correlator.\n";
        bit_image = true;
      } else if (vm.count("nonbitimage")) {
        std::cout << "Forcing the use of the non bit-image optimized correlator.\n";
        bit_image = false;
      } 

      if (bit_image) {
        ImageView<PixelGray<uint8> > left_bitimage = channel_cast<uint8>(left);
        ImageView<PixelGray<uint8> > right_bitimage = channel_cast<uint8>(right);
        vw::Timer corr_timer("Correlation Time");
        disparity_map = correlator( left_bitimage, right_bitimage, bit_image );
      } else {
        vw::Timer corr_timer("Correlation Time");
        disparity_map = correlator( left, right, bit_image );
      }
    }

    remove_outliers( disparity_map, 5, 5, 3, 0.6 );

    // Write disparity debug images
    BBox2 disp_range = disparity::get_disparity_range(disparity_map);
    write_image( "x_disparity.png", normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x() )));
    write_image( "y_disparity.png", normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y() )));
  }
  catch( vw::Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}

