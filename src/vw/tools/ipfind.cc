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


/// \file ipfind.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///

#include <vw/Core/System.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralDescriptor.h>
#include <vw/InterestPoint/IntegralDetector.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>
#include <vw/InterestPoint/InterestData.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;
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
  float inc_amt = 1/norm_2(delta);
  for ( float r=0; r<1.0; r+=inc_amt ) {
    int i = (int)(0.5 + start.x() + r*float(delta.x()) );
    int j = (int)(0.5 + start.y() + r*float(delta.y()) );
    image.impl()(i,j) = value;
  }
}

static void write_debug_image( std::string const& out_file_name,
                               std::string const& input_file_name,
                               InterestPointList const& ip ) {
  vw_out() << "Writing debug image: " << out_file_name << "\n";

  vw_out(InfoMessage,"interest_point") << "\t > Gathering statistics:\n";
  float min = 1e30, max = -1e30;
  for ( InterestPointList::const_iterator point = ip.begin(); point != ip.end(); ++point ) {
    if ( point->interest > max )
      max = point->interest;
    if ( point->interest < min )
      min = point->interest;
  }
  float diff = max - min;

  vw_out(InfoMessage,"interest_point") << "\t > Drawing raster:\n";
  ImageView<PixelRGB<uint8> > oimage;
  boost::scoped_ptr<DiskImageResource> irsrc( DiskImageResource::open(input_file_name) );
  if ( irsrc->has_nodata_read() ) {
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(apply_mask(normalize(create_mask(DiskImageView<PixelGray<float> >( *irsrc ),
                                                                                   irsrc->nodata_read()))) * 0.5 );
  } else {
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(normalize(DiskImageView<PixelGray<float> >( *irsrc )) * 0.5);
  }

  for ( InterestPointList::const_iterator point = ip.begin(); point != ip.end(); ++point ) {
    float norm_i = (point->interest - min)/diff;
    PixelRGB<uint8> color(0,0,0);
    if ( norm_i < .5 ) {
      // Form of red
      color.r() = 255;
      color.g() = (uint8)(2*norm_i*255);
    } else {
      // Form of green
      color.g() = 255;
      color.r() = 255 - (uint8)(2*(norm_i-.5)*255);
    }

    // Marking point w/ Dot
    oimage(point->ix,point->iy) = color;

    // Circling point based on the scale
    for (float a = 0; a < 6; a+=.392 ) {
      float a_d = a + .392;
      Vector2i start( int(2*point->scale*cos(a  )+point->x),
                      int(2*point->scale*sin(a  )+point->y) );
      Vector2i end(   int(2*point->scale*cos(a_d)+point->x),
                      int(2*point->scale*sin(a_d)+point->y) );
      draw_line( oimage, color, start, end );
    }
  }

  boost::scoped_ptr<DiskImageResource> rsrc( DiskImageResource::create(out_file_name,
                                                                       oimage.format() ) );
  vw_out(InfoMessage,"interest_point") << "\t > Writing out image:\n";
  block_write_image( *rsrc, oimage,
                     TerminalProgressCallback( "tools.ipfind","\t : ") );
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  std::string interest_operator, descriptor_generator;
  float  ip_gain;
  uint32 max_points;
  int    tile_size, num_threads;
  ImageView<double> integral;
  bool   no_orientation;

  const float IDEAL_LOG_THRESHOLD    = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("num-threads", po::value(&num_threads)->default_value(0), 
          "Set the number of threads for interest point detection.  Setting the num_threads to zero causes ipfind to use the visionworkbench default number of threads.")
    ("tile-size,t", po::value(&tile_size), 
          "Specify the tile size for processing interest points. (Useful when working with large images).")
    ("lowe,l",        "Save the interest points in an ASCII data format that is compatible with the Lowe-SIFT toolchain.")
    ("normalize",     "Normalize the input, use for images that have non standard values such as ISIS cube files.")
    ("debug-image,d", "Write out debug images.")

    // Interest point detector options
    ("interest-operator", po::value(&interest_operator)->default_value("IAGD"), 
          "Choose an interest point metric from [IAGD (ASP default), LoG, Harris, OBALoG, Brisk, Orb]")
    ("gain,g", po::value(&ip_gain)->default_value(1.0), 
          "Increasing this number will increase the gain at which interest points are detected.")
    ("max-points", po::value(&max_points)->default_value(0), 
          "Set the maximum number of interest points you want returned.  The most \"interesting\" points are selected.")
    ("single-scale", 
          "Turn off scale-invariant interest point detection.  This option only searches for interest points in the first octave of the scale space.")
    ("no-orientation", po::bool_switch(&no_orientation), "Shutoff rotational invariance")

    // Descriptor generator options
    ("descriptor-generator", po::value(&descriptor_generator)->default_value("sgrad"), 
          "Choose a descriptor generator from [patch, pca, sgrad (ASP default), sgrad2]");

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
  } catch (const po::error& e) {
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

  if ( num_threads > 0 )
    vw_settings().set_default_num_threads(num_threads);
  if ( vm.count("tile-size"))
    vw_settings().set_default_tile_size(tile_size);

  // Checking strings
  boost::to_lower( interest_operator );
  boost::to_lower( descriptor_generator );
  // Determine if interest_operator is legitimate
  if ( !( interest_operator == "iagd"   ||
          interest_operator == "harris" ||
          interest_operator == "log"    ||
          interest_operator == "obalog" || 
          interest_operator == "brisk"  ||
          interest_operator == "orb") ) {
    vw_out() << "Unknown interest operator: " << interest_operator
             << ". Options are : [ IAGD, Harris, LoG, OBALoG, Brisk, Orb ]\n";
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
  for (size_t i = 0; i < input_file_names.size(); ++i) {

    vw_out() << "Finding interest points in \"" << input_file_names[i] << "\".\n";
    std::string file_prefix = fs::path(input_file_names[i]).replace_extension().string();
    boost::scoped_ptr<DiskImageResource> image_rsrc( DiskImageResource::open( input_file_names[i] ) );
    DiskImageView<PixelGray<float> > raw_image( *image_rsrc );
    ImageViewRef<PixelGray<float> >  image = raw_image;
    if ( vm.count("normalize") && image_rsrc->has_nodata_read() )
      image = apply_mask(normalize(create_mask(raw_image,image_rsrc->nodata_read())));
    else if ( vm.count("normalize") )
      image = normalize(raw_image);

    // Potentially mask image on a no data value
    if ( image_rsrc->has_nodata_read() )
      vw_out(DebugMessage,"interest_point") << "Image has a nodata value: "
                                            << image_rsrc->nodata_read() << "\n";

    // The max points sent to IP Detector class is applied to each
    // tile of an image. In order to curb memory use we'll set the max
    // size for each tile smaller (proportional to the number of tiles).
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
        ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator, tile_max_points);
        ip = detect_interest_points(image, detector);
      } else {
        InterestPointDetector<HarrisInterestOperator> detector(interest_operator, tile_max_points);
        ip = detect_interest_points(image, detector);
      }
    } else if ( interest_operator == "log") {
      // Use a scale-space Laplacian of Gaussian feature detector. The
      // associated threshold is abs(interest) > interest_threshold.
      // LoG threshold is inversely proportional to gain..
      LogInterestOperator interest_operator(IDEAL_LOG_THRESHOLD/ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator, tile_max_points);
        ip = detect_interest_points(image, detector);
      } else {
        InterestPointDetector<LogInterestOperator> detector(interest_operator, tile_max_points);
        ip = detect_interest_points(image, detector);
      }
    } else if ( interest_operator == "obalog") {
      // OBALoG threshold is inversely proportional to gain ..
      OBALoGInterestOperator interest_operator(IDEAL_OBALOG_THRESHOLD/ip_gain);
      IntegralInterestPointDetector<OBALoGInterestOperator> detector( interest_operator, tile_max_points );
      ip = detect_interest_points(image, detector);
    } else if ( interest_operator == "iagd") {
      // This is the default ASP implementation
      IntegralAutoGainDetector detector( tile_max_points );
      ip = detect_interest_points(image, detector);
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    } else if ( interest_operator == "brisk") {
      //
      OpenCvInterestPointDetector detector(OPENCV_IP_DETECTOR_TYPE_BRISK);
      ip = detect_interest_points(image, detector);
    } else if ( interest_operator == "orb") {
      //
      OpenCvInterestPointDetector detector(OPENCV_IP_DETECTOR_TYPE_ORB);
      ip = detect_interest_points(image, detector);
    }
#else // End OpenCV section
    } else {
      vw_out() << "Error: Cannot use ORB or BRISK if not compiled with OpenCV!" << std::endl;
      exit(0); 
    }
#endif

    // Removing Interest Points on nodata or within 1/px
    if ( image_rsrc->has_nodata_read() ) {
      ImageViewRef<PixelMask<PixelGray<float> > > image_mask =
        create_mask(raw_image,image_rsrc->nodata_read());
      BBox2i bound = bounding_box( image_mask );
      bound.contract(1);
      int before_size = ip.size();
      for ( InterestPointList::iterator point = ip.begin(); point != ip.end(); ++point ) {

        // To Avoid out of index issues
        if ( !bound.contains( Vector2i(point->ix, point->iy ))) {
          point = ip.erase(point);
          point--;
          continue;
        }

        if ( !image_mask(point->ix,  point->iy  ).valid() ||
             !image_mask(point->ix+1,point->iy+1).valid() ||
             !image_mask(point->ix+1,point->iy  ).valid() ||
             !image_mask(point->ix+1,point->iy-1).valid() ||
             !image_mask(point->ix,  point->iy+1).valid() ||
             !image_mask(point->ix,  point->iy-1).valid() ||
             !image_mask(point->ix-1,point->iy+1).valid() ||
             !image_mask(point->ix-1,point->iy  ).valid() ||
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
    if ( max_points > 0  && ip.size() > max_points ) {
      ip.sort();
      ip.resize(max_points);
      vw_out() << "\t Culled to " << ip.size() << " points.\n";
    }

    // Delete the orientation value if requested
    if ( no_orientation ) {
      BOOST_FOREACH( InterestPoint& i, ip ) {
        i.orientation = 0;
      }
    }

    // Generate descriptors for interest points.
    vw_out(InfoMessage) << "\tRunning " << descriptor_generator << " descriptor generator.\n";
    if (descriptor_generator == "patch") {
      PatchDescriptorGenerator descriptor;
      describe_interest_points( image, descriptor, ip );
    } else if (descriptor_generator == "pca") {
      PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
      describe_interest_points( image, descriptor, ip );
    } else if (descriptor_generator == "sgrad") {
      SGradDescriptorGenerator descriptor;
      describe_interest_points( image, descriptor, ip );
    } else if (descriptor_generator == "sgrad2") {
      SGrad2DescriptorGenerator descriptor;
      describe_interest_points( image, descriptor, ip );
    }

    // If ASCII output was requested, write it out.  Otherwise stick with binary output.
    if (vm.count("lowe"))
      write_lowe_ascii_ip_file(file_prefix + ".key", ip);
    else
      write_binary_ip_file(file_prefix + ".vwip", ip);

    // Write Debug image
    if (vm.count("debug-image"))
      write_debug_image( file_prefix + "_debug.jpg",
                         input_file_names[i],
                         ip );
  }

  return 0;
}
