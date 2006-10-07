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

    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("left", po::value<std::string>(&left_file_name), "Explicitly specify the \"left\" input file")
      ("right", po::value<std::string>(&right_file_name), "Explicitly specify the \"right\" input file")
      ("slog", po::value<float>(&slog)->default_value(0.0), "Apply SLOG filter with the given sigma, or 0 to disable")
      ("log", po::value<float>(&log)->default_value(0.0), "Apply LOG filter with the given sigma, or 0 to disable")
      ("xoffset", po::value<int>(&xoffset)->default_value(0), "Overall horizontal offset between images")
      ("yoffset", po::value<int>(&yoffset)->default_value(0), "Overall vertical offset between images")
      ("xrange", po::value<int>(&xrange)->default_value(5), "Allowed range of horizontal disparity")
      ("yrange", po::value<int>(&yrange)->default_value(5), "Allowed range of vertical disparity")
      ("xkernel", po::value<int>(&xkernel)->default_value(15), "Horizontal correlation kernel size")
      ("ykernel", po::value<int>(&ykernel)->default_value(15), "Vertical correlation kernel size")
      ("lrthresh", po::value<int>(&lrthresh)->default_value(2), "Left/right correspondence threshold")
      ("hsubpix", "Enable horizontal sub-pixel correlation")
      ("vsubpix", "Enable vertical sub-pixel correlation")
      ("reference", "Use the slower, simpler reference correlator")
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
    if( log > 0 ) {
      std::cout << "Applying LOG filter..." << std::endl;
      left = laplacian_filter( gaussian_filter( left, log ) );
      right = laplacian_filter( gaussian_filter( right, log ) );
    }

    if( slog > 0.0 ) {
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
    } else {
      vw::stereo::SubpixelCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                 yoffset-yrange, yoffset+yrange,
                                                 xkernel, ykernel,
                                                 true, lrthresh,
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
        disparity_map = correlator( left_bitimage, right_bitimage, bit_image );
      } else {
        disparity_map = correlator( left, right, bit_image );
      }
    }

    double min_horz_disp, max_horz_disp, min_vert_disp, max_vert_disp;
    get_disparity_range( disparity_map, 
                         min_horz_disp, max_horz_disp, 
                         min_vert_disp, max_vert_disp,
                         true);

    remove_outliers( disparity_map, 5, 5, 60, 3 );

    ImageView<PixelGray<float> > x_disparity, y_disparity;
    disparity_debug_images( disparity_map, x_disparity, y_disparity );
    write_image( "x_disparity.png", normalize(x_disparity) );
    write_image( "y_disparity.png", normalize(y_disparity) );
  }
  catch( vw::Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}
