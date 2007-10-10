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
#include <vw/Stereo/PyramidCorrelator.h>

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
      ("hsubpix", "Enable horizontal sub-pixel correlation")
      ("vsubpix", "Enable vertical sub-pixel correlation")
      ("reference", "Use the slower, simpler reference correlator")
      ("pyramid", "Use the pyramid based correlator")
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

    DiskImageView<PixelGray<float> > left_disk_image(left_file_name );
    DiskImageView<PixelGray<float> > right_disk_image(right_file_name );
    int cols = std::max(left_disk_image.cols(),right_disk_image.cols());
    int rows = std::max(left_disk_image.rows(),right_disk_image.rows());
    ImageViewRef<PixelGray<float> > left = edge_extend(left_disk_image,0,0,cols,rows);
    ImageViewRef<PixelGray<float> > right = edge_extend(right_disk_image,0,0,cols,rows);

    ImageView<PixelDisparity<float> > disparity_map;
    if (vm.count("reference")>0) {
      vw::stereo::ReferenceCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                  yoffset-yrange, yoffset+yrange,
                                                  xkernel, ykernel,
                                                  true, lrthresh,
                                                  (vm.count("hsubpix")>0),
                                                  (vm.count("vsubpix")>0) );
      if (log > 0) 
        disparity_map = correlator( left, right, stereo::LogStereoPreprocessingFilter(log));
      else 
        disparity_map = correlator( left, right, stereo::SlogStereoPreprocessingFilter(slog));
    } else if (vm.count("pyramid")>0) {
      vw::stereo::PyramidCorrelator correlator( BBox2(Vector2(xoffset-xrange, yoffset-yrange),
                                                      Vector2(xoffset+xrange, yoffset+yrange)),
                                                Vector2i(xkernel, ykernel),
                                                slog, lrthresh,
                                                corrscore_thresh,
                                                (vm.count("hsubpix")>0),
                                                (vm.count("vsubpix")>0) );
      correlator.set_debug_mode("debug");
      if (1) {
        vw::Timer corr_timer("Correlation Time");
        if (log > 0) 
          disparity_map = correlator( left, right, stereo::LogStereoPreprocessingFilter(log));
         else 
           disparity_map = correlator( left, right, stereo::SlogStereoPreprocessingFilter(slog));
      }
    } else {
      vw::stereo::OptimizedCorrelator correlator( xoffset-xrange, xoffset+xrange,
                                                  yoffset-yrange, yoffset+yrange,
                                                  xkernel, ykernel,
                                                  true, lrthresh, corrscore_thresh,
                                                  (vm.count("hsubpix")>0),
                                                  (vm.count("vsubpix")>0) );
      if (1) {
        vw::Timer corr_timer("Correlation Time");
        if (log > 0) 
          disparity_map = correlator( left, right, stereo::LogStereoPreprocessingFilter(log));
        else 
          disparity_map = correlator( left, right, stereo::SlogStereoPreprocessingFilter(slog));
      }
    }

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

