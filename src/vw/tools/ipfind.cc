// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ipfind.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///
#include <vw/Core.h>
#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Math/BBox.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

template <class ImageT, class ValueT>
void draw_line( ImageViewBase<ImageT>& image,
                ValueT const& value,
                Vector2i const& start,
                Vector2i const& end ) {

  BBox2i bound = bounding_box(image.impl());
  if ( !bound.contains( start ) ||
       !bound.contains( end ) )
    return;
  Vector2i delta = end - start;
  for ( float r=0; r<1.0; r+=1/norm_2(delta) ) {
    int i = (int)(0.5 + start.x() + r*float(delta.x()) );
    int j = (int)(0.5 + start.y() + r*float(delta.y()) );
    image.impl()(i,j) = value;
  }
}

static void write_debug_image( std::string out_file_name,
                               std::string input_file_name,
                               InterestPointList const& ip ) {
  vw_out() << "Writing debug image: " << out_file_name << "\n";
  DiskImageView<PixelGray<uint8> > image( input_file_name );

  vw_out(InfoMessage,"interest_point") << "\t > Gathering statistics:\n";
  float min = 1e30, max = -1e30;
  for ( InterestPointList::const_iterator point = ip.begin();
        point != ip.end(); ++point ) {
    if ( point->interest > max )
      max = point->interest;
    if ( point->interest < min )
      min = point->interest;
  }
  float diff = max - min;

  vw_out(InfoMessage,"interest_point") << "\t > Drawing raster:\n";
  ImageView<PixelRGB<uint8> > oimage;
  oimage = pixel_cast<PixelRGB<uint8> >(image*0.5);
  for ( InterestPointList::const_iterator point = ip.begin();
        point != ip.end(); ++point ) {
    float norm_i = (point->interest - min)/diff;
    PixelRGB<uint8> color(0,0,0);
    if ( norm_i < .5 ) {
      // Form of red
      color.r() = 255;
      color.g() = (unsigned char)(2*norm_i*255);
    } else {
      // Form of green
      color.g() = 255;
      color.r() = 255 - (unsigned char)(2*(norm_i-.5)*255);
    }

    // Marking point w/ Dot
    oimage(point->ix,point->iy) = color;

    // Circling point
    for (float a = 0; a < 6; a+=.392 ) {
      float a_d = a + .392;
      Vector2i start( int(2*point->scale*cos(a)+point->x),
                      int(2*point->scale*sin(a)+point->y) );
      Vector2i end( int(2*point->scale*cos(a_d)+point->x),
                    int(2*point->scale*sin(a_d)+point->y) );
      draw_line( oimage, color, start, end );
    }
  }

  DiskImageResource *rsrc = DiskImageResource::create(out_file_name,
                                                      oimage.format() );
  vw_out(InfoMessage,"interest_point") << "\t > Writing out image:\n";
  block_write_image( *rsrc, oimage,
                     TerminalProgressCallback( "tools.ipfind","\t : ") );

}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  std::string interest_operator, descriptor_generator;
  float ip_gain;
  uint32 max_points;
  int tile_size, num_threads;
  ImageView<double> integral;

  const float IDEAL_LOG_THRESHOLD = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("num-threads", po::value(&num_threads)->default_value(0), "Set the number of threads for interest point detection.  Setting the num_threads to zero causes ipfind to use the visionworkbench default number of threads.")
    ("tile-size,t", po::value(&tile_size), "Specify the tile size for processing interest points. (Useful when working with large images). VW usually picks 1024 px.")
    ("lowe,l", "Save the interest points in an ASCII data format that is compatible with the Lowe-SIFT toolchain.")
    ("debug-image,d", "Write out debug images.")

    // Interest point detector options
    ("interest-operator", po::value(&interest_operator)->default_value("OBALoG"), "Choose an interest point metric from [LoG, Harris, OBALoG]")
    ("gain,g", po::value(&ip_gain)->default_value(1.0), "Increasing this number will increase the gain at which interest points are detected.")
    ("max-points", po::value(&max_points)->default_value(0), "Set the maximum number of interest points you want returned.  The most \"interesting\" points are selected.")
    ("single-scale", "Turn off scale-invariant interest point detection.  This option only searches for interest points in the first octave of the scale space.")

    // Descriptor generator options
    ("descriptor-generator", po::value(&descriptor_generator)->default_value("sgrad"), "Choose a descriptor generator from [patch,pca,sgrad,sgrad2]");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filenames>..." << std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e ) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    vw_out() << usage.str();
    return 1;
  }

  if( input_file_names.size() < 1 ) {
    vw_out() << "Error: Must specify at least one input file!" << std::endl << std::endl;
    vw_out() << usage.str();
    return 1;
  }

  if ( vm.count("num-threads"))
    vw_settings().set_default_num_threads(num_threads);
  if ( vm.count("tile-size"))
    vw_settings().set_default_tile_size(tile_size);

  // Checking strings
  boost::to_lower( interest_operator );
  boost::to_lower( descriptor_generator );
  // Determine if interest_operator is legitimate
  if ( !( interest_operator == "harris" ||
          interest_operator == "log" ||
          interest_operator == "obalog" ) ) {
    vw_out() << "Unknown interest operator: " << interest_operator
              << ". Options are : [ Harris, LoG, OBALoG ]\n";
    exit(0);
  }
  // Determine if descriptor_generator is legitimate
  if ( !( descriptor_generator == "patch" ||
          descriptor_generator == "pca"   ||
          descriptor_generator == "sgrad" ||
          descriptor_generator == "sgrad2" ) ) {
    vw_out() << "Unkown descriptor generator: " << descriptor_generator
              << ". Options are : [ Patch, PCA, SGrad, SGrad2 ]\n";
    exit(0);
  }

  // Iterate over the input files and find interest points in each.
  for (unsigned i = 0; i < input_file_names.size(); ++i) {

    vw_out() << "Finding interest points in \"" << input_file_names[i] << "\".\n";
    std::string file_prefix = fs::path(input_file_names[i]).replace_extension().string();
    DiskImageResource *image_rsrc = DiskImageResource::open( input_file_names[i] );
    DiskImageView<PixelGray<float> > image(image_rsrc);

    // Potentially mask image on a no data value
    if ( image_rsrc->has_nodata_read() )
      vw_out(DebugMessage,"interest_point") << "Image has a nodata value: "
                                            << image_rsrc->nodata_read() << "\n";

    // The max points sent to IP Detector class is applied to each
    // tile of an image. In order to curb memory use we'll set the max
    // size for each tile smaller (proportional to the number of
    // tiles).
    int number_tiles = (image.cols()/vw_settings().default_tile_size()+1) *
      (image.rows()/vw_settings().default_tile_size()+1);
    uint32 tile_max_points = uint32(float(max_points)/float(number_tiles))*2; // A little over shoot
                                                                     // incase the tile is empty
    if ( max_points == 0 ) tile_max_points = 0; // No culling
    else if ( tile_max_points < 50 ) tile_max_points = 50;

    // Detecting Interest Points
    InterestPointList ip;
    if ( interest_operator == "harris" ) {
      // Harris threshold is inversely proportional to gain.
      HarrisInterestOperator interest_operator(IDEAL_HARRIS_THRESHOLD/ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator,
                                                                     tile_max_points);
        ip = detect_interest_points(image, detector);
      } else {
        InterestPointDetector<HarrisInterestOperator> detector(interest_operator,
                                                               tile_max_points);
        ip = detect_interest_points(image, detector);
      }
    } else if ( interest_operator == "log") {
      // Use a scale-space Laplacian of Gaussian feature detector. The
      // associated threshold is abs(interest) > interest_threshold.
      // LoG threshold is inversely proportional to gain..
      LogInterestOperator interest_operator(IDEAL_LOG_THRESHOLD/ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator,
                                                                  tile_max_points);
        ip = detect_interest_points(image, detector);
      } else {
        InterestPointDetector<LogInterestOperator> detector(interest_operator,
                                                            tile_max_points);
        ip = detect_interest_points(image, detector);
      }
    } else if ( interest_operator == "obalog") {
      // OBALoG threshold is inversely proportional to gain ..
      OBALoGInterestOperator interest_operator(IDEAL_OBALOG_THRESHOLD/ip_gain);
      IntegralInterestPointDetector<OBALoGInterestOperator> detector( interest_operator,
                                                                      tile_max_points );
      ip = detect_interest_points(image, detector);
    }

    // Removing Interest Points on nodata or within 1/px
    if ( image_rsrc->has_nodata_read() ) {
      float nodata_value = image_rsrc->nodata_read();
      ImageViewRef<PixelMask<PixelGray<float> > > image_mask = create_mask(image,nodata_value);
      BBox2i bound = bounding_box( image_mask );
      bound.contract(1);
      int before_size = ip.size();
      for ( InterestPointList::iterator point = ip.begin();
            point != ip.end(); ++point ) {

        // To Avoid out of index issues
        if ( !bound.contains( Vector2i(point->ix,
                                       point->iy ))) {
          point = ip.erase(point);
          point--;
          continue;
        }

        if ( !image_mask(point->ix,point->iy).valid() ||
             !image_mask(point->ix+1,point->iy+1).valid() ||
             !image_mask(point->ix+1,point->iy).valid() ||
             !image_mask(point->ix+1,point->iy-1).valid() ||
             !image_mask(point->ix,point->iy+1).valid() ||
             !image_mask(point->ix,point->iy-1).valid() ||
             !image_mask(point->ix-1,point->iy+1).valid() ||
             !image_mask(point->ix-1,point->iy).valid() ||
             !image_mask(point->ix-1,point->iy-1).valid() ) {
          point = ip.erase(point);
          point--;
          continue;
        }
      }
      vw_out(InfoMessage,"interest_point") << "Removed " << before_size-ip.size() << " points close to nodata.\n";
    }

    vw_out() << "\t Found " << ip.size() << " points.\n";

    // Additional Culling for the entire image
    ip.sort();
    if ( (max_points > 0) && (ip.size() > max_points) ) {
      ip.resize(max_points);
      vw_out() << "\t Culled to " << ip.size() << " points.\n";
    }

    // Generate descriptors for interest points.
    vw_out(InfoMessage) << "\tRunning " << descriptor_generator << " descriptor generator.\n";
    if (descriptor_generator == "patch") {
      PatchDescriptorGenerator descriptor;
      descriptor(image, ip);
    } else if (descriptor_generator == "pca") {
      PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
      descriptor(image, ip);
    } else if (descriptor_generator == "sgrad") {
      SGradDescriptorGenerator descriptor;
      descriptor(image, ip);
    } else if (descriptor_generator == "sgrad2") {
      SGrad2DescriptorGenerator descriptor;
      descriptor(image, ip);
    }

    // If ASCII output was requested, write it out.  Otherwise stick
    // with binary output.
    if (vm.count("lowe"))
      write_lowe_ascii_ip_file(file_prefix + ".key", ip);
    else
      write_binary_ip_file(file_prefix + ".vwip", ip);

    // Write Debug image
    if (vm.count("debug-image")) {
      std::string output_file_name =
        file_prefix + "_debug.jpg";
      write_debug_image( output_file_name,
                         input_file_names[i],
                         ip );
    }
  }
}
