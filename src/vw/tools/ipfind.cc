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

// TODO: Move this to the IP library, then remove from ASP.
/// Remove all pixels within a radius of the border or nodata.
// Note: A nodata pixel is one for which pixel <= nodata.
template <class ImageT>
void remove_ip_near_nodata( vw::ImageViewBase<ImageT> const& image,   double nodata,
                            vw::ip::InterestPointList      & ip_list, int    radius ){

  size_t prior_ip = ip_list.size();
  
  typedef ImageView<typename ImageT::pixel_type> CropImageT;
  const int width = 2*radius+1;
  CropImageT subsection(width,width);   
  
  // Get shrunk bounding box
  BBox2i bound = bounding_box( image.impl() );
  bound.contract(radius); 
  // Loop through all the points
  for ( ip::InterestPointList::iterator ip = ip_list.begin(); ip != ip_list.end(); ++ip ) {
    
    // Remove the point if it was too close to the image borders
    if ( !bound.contains( Vector2i(ip->ix,ip->iy) ) ) {
      ip = ip_list.erase(ip);
      ip--;
      continue;
    }

    // Remove the point if any invalid pixels are too close to it
    subsection = crop( image.impl(), ip->ix-radius, ip->iy-radius, width, width );
    for ( typename CropImageT::iterator pixel = subsection.begin();
          pixel != subsection.end(); pixel++ ) {
      if (*pixel <= nodata) {
        ip = ip_list.erase(ip);
        ip--;
        break;
      }
    }
  }
  VW_OUT( DebugMessage, "asp" ) << "Removed " << prior_ip - ip_list.size()
			  << " interest points due to their proximity to nodata values."
			  << std::endl << "Nodata value used " << nodata << std::endl;
} // End function remove_ip_near_nodata

// TODO: This should be merged with the debug functions in ASP

/// Draw a helpful image showing where IPs are detected in an image.
/// - If reduce_scale is set the IP circle markers will be smaller.
template <typename ImageT>
void write_debug_image( std::string const& out_file_name,
                        vw::ImageViewBase<ImageT> const& image,
                        InterestPointList const& ip,
                        bool has_nodata, double nodata,
                        bool force_full_res=false,
                        bool reduce_scale=false) {

  // Scale the images to keep the size down below 1500x1500 so they draw quickly.
  float sub_scale  = sqrt(1500.0 * 1500.0 / float(image.impl().cols() * image.impl().rows()));
	sub_scale /= 2;
  if ((sub_scale > 1) || force_full_res)
    sub_scale = 1;

  vw_out() << "Writing debug image: " << out_file_name << " with downsample: "<< sub_scale << "\n";
  
  // Get the min and max point intensity
  vw_out(InfoMessage,"interest_point") << "\t > Gathering statistics:\n";
  float min = 1e30, max = -1e30;
  for ( InterestPointList::const_iterator point = ip.begin(); point != ip.end(); ++point ) {
    if ( point->interest > max )
      max = point->interest;
    if ( point->interest < min )
      min = point->interest;
  }
  float diff = max - min;

  // Create a muted version of the input image
  ImageView<PixelRGB<uint8> > oimage;
  const double IMAGE_FADE_PERCENTAGE = 0.7;
  if ( has_nodata ) {
    ImageView<PixelMask<float> > gray = resample(create_mask(image, nodata), sub_scale);
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(apply_mask(normalize(gray) * IMAGE_FADE_PERCENTAGE));
  } else {
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(normalize(resample(image, sub_scale)) * IMAGE_FADE_PERCENTAGE);
  }

  // Draw each of the interest points on to the muted input image.
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
    int ix = point->ix*sub_scale;
    int iy = point->iy*sub_scale;
    int x  = point->x*sub_scale;
    int y  = point->y*sub_scale;
    oimage(ix,iy) = color;
    
    double scale = point->scale;
    if (reduce_scale)
      scale = 0.12*scale;

    // Circling point based on the scale
    for (float a = 0; a < 6; a+=.392 ) {
      float a_d = a + .392; // ?
      Vector2i start( int(2*scale*cos(a  )+x),
                      int(2*scale*sin(a  )+y) );
      Vector2i end(   int(2*scale*cos(a_d)+x),
                      int(2*scale*sin(a_d)+y) );
      draw_line( oimage, color, start, end );
    }
  } // End loop through interest points

  boost::scoped_ptr<DiskImageResource> rsrc( DiskImageResource::create(out_file_name,
                                                                       oimage.format() ) );
  vw_out(InfoMessage,"interest_point") << "\t > Writing out image:\n";
  write_image( *rsrc, oimage,
                     TerminalProgressCallback( "tools.ipfind","\t : ") );
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  std::string output_folder, interest_operator, descriptor_generator;
  float  ip_gain;
  uint32 max_points;
  int    tile_size, num_threads, nodata_radius, print_num_ip, debug_image;
  ImageView<double> integral;
  bool   no_orientation;
  bool   opencv_normalize = false;

  const float IDEAL_LOG_THRESHOLD    = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")
    ("output-folder",  po::value(&output_folder)->default_value(""), 
                      "Write output files to this location.")
    ("num-threads",          po::value(&num_threads)->default_value(0), 
                      "Set the number of threads for interest point detection.  Setting the num_threads to zero causes ipfind to use the visionworkbench default number of threads.")
    ("tile-size,t",          po::value(&tile_size), 
                      "Specify the tile size for processing interest points. (Useful when working with large images).")
    ("lowe",          "Save the interest points in an ASCII data format that is compatible with the Lowe-SIFT toolchain.")
    ("normalize",     "Normalize the input, use for images that have non standard values such as ISIS cube files.")
    ("opencv-normalize",     "Apply per-tile openCV normalization.")
    ("nodata-radius,", po::value(&nodata_radius)->default_value(1), 
                      "Don't detect IP within this many pixels of image borders or nodata.")
    ("debug-image,d", po::value(&debug_image)->default_value(0), 
                      "Write out debug images (1), write out full res debug images (2).")
    ("print-ip",             po::value(&print_num_ip)->default_value(0), 
                      "Print information for this many detected IP")

    // Interest point detector options
    ("interest-operator",    po::value(&interest_operator)->default_value("sift"), 
                      "Choose an interest point metric from [IAGD (ASP default), LoG, Harris, OBALoG, Orb, Sift]")
    ("gain,g",               po::value(&ip_gain)->default_value(1.0), 
                      "Increasing this number will increase the gain at which interest points are detected.")
    ("max-points",           po::value(&max_points)->default_value(250), 
                      "Set the maximum number of interest points you want returned per tile.")
    ("single-scale", 
                      "Turn off scale-invariant interest point detection.  This option only searches for interest points in the first octave of the scale space. Harris and LoG only.")
    ("no-orientation",       po::bool_switch(&no_orientation), "Shutoff rotational invariance")

    // Descriptor generator options
    ("descriptor-generator", po::value(&descriptor_generator)->default_value("sift"), 
                      "Choose a descriptor generator from [patch, pca, sgrad, sift, sgrad2, orb]");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  std::ostringstream usage;
  usage << "Usage: ipfind [options] <filenames>...\n\n";
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

  if( nodata_radius < 1 ) {
    vw_out() << "Error: nodata-radius must be >= 1!" << std::endl << std::endl;
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

  // The OpenCV detector types are handled differently than the custom detector types.  
  const bool detector_is_opencv = ( (interest_operator == "brisk") ||
                                    (interest_operator == "surf" ) ||
                                    (interest_operator == "orb"  ) ||
                                    (interest_operator == "sift" )   );
  
  // Determine if interest_operator is legitimate
  if ( !( interest_operator == "iagd"   ||
          interest_operator == "harris" ||
          interest_operator == "log"    ||
          interest_operator == "obalog" || 
          //interest_operator == "brisk"  ||
          interest_operator == "orb"    ||
          //interest_operator == "surf"   ||
          interest_operator == "sift"     ) ) {
    vw_out() << "Unknown interest operator: " << interest_operator
             << ". Options are : [ IAGD, Harris, LoG, OBALoG, Orb, Sift ]\n";
    exit(0);
  }
  // Determine if descriptor_generator is legitimate
  if ( !( descriptor_generator == "patch"  ||
          descriptor_generator == "pca"    ||
          descriptor_generator == "sgrad"  ||
          descriptor_generator == "sgrad2" ||
          //descriptor_generator == "brisk"  ||
          descriptor_generator == "orb"    ||
          //descriptor_generator == "surf"   ||
          descriptor_generator == "sift"    ) ) {
    vw_out() << "Unkown descriptor generator: " << descriptor_generator
             << ". Options are : [ Patch, PCA, SGrad, SGrad2, Orb, Sift ]\n";
    exit(0);
  }

  if ( vm.count("opencv-normalize"))
    opencv_normalize = true;

  // Iterate over the input files and find interest points in each.
  for (size_t i = 0; i < input_file_names.size(); ++i) {

    vw_out() << "Finding interest points in \"" << input_file_names[i] << "\".\n";
    std::string file_prefix = boost::filesystem::path(input_file_names[i]).replace_extension().string();
    if (output_folder != "") {
      // Write the output files to the specified output folder.
      boost::filesystem::path of(output_folder);
      of /= boost::filesystem::path(file_prefix).filename();
      file_prefix = of.string();
    }

    // Get reference to the file and convert to floating point
    boost::shared_ptr<vw::DiskImageResource> image_rsrc( vw::DiskImageResourcePtr( input_file_names[i] ) );
    DiskImageView<PixelGray<float> > raw_image( *image_rsrc );
    ImageViewRef<PixelGray<float> >  image;
    ImageViewRef<PixelMask<PixelGray<float> > >  masked_image;

    // Load nodata value if present.
    double nodata = 0;
    bool   has_nodata = image_rsrc->has_nodata_read();
    if (has_nodata)
      nodata = image_rsrc->nodata_read();

    if ( vm.count("normalize") && has_nodata ) {
      masked_image = normalize(create_mask(raw_image, nodata));
      image        = apply_mask(masked_image); // A backup option for detectors which don't handle nodata.
    } else {
      if ( vm.count("normalize") )
        image = normalize(raw_image);
      else
        image = raw_image;
    }

    // Potentially mask image on a no data value
    if (has_nodata)
      vw_out(DebugMessage,"interest_point") << "Image has a nodata value: " << nodata << "\n";

    const bool describeInDetect = true;

    // Detecting Interest Points
    // - Note that only the OpenCV detectors handle masks.
    InterestPointList ip;
    if ( interest_operator == "harris" ) {
      // Harris threshold is inversely proportional to gain.
      HarrisInterestOperator interest_operator(IDEAL_HARRIS_THRESHOLD/ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, max_points);
      } else {
        InterestPointDetector<HarrisInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, max_points);
      }
    } else if ( interest_operator == "log") {
      // Use a scale-space Laplacian of Gaussian feature detector. The
      // associated threshold is abs(interest) > interest_threshold.
      // LoG threshold is inversely proportional to gain..
      LogInterestOperator interest_operator(IDEAL_LOG_THRESHOLD/ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, max_points);
      } else {
        InterestPointDetector<LogInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, max_points);
      }
    } else if ( interest_operator == "obalog") {
      // OBALoG threshold is inversely proportional to gain ..
      OBALoGInterestOperator interest_operator(IDEAL_OBALOG_THRESHOLD/ip_gain);
      IntegralInterestPointDetector<OBALoGInterestOperator> detector( interest_operator, max_points );
      ip = detect_interest_points(image, detector, max_points);
    } else if ( interest_operator == "iagd") {
      // This is the default ASP implementation
      IntegralAutoGainDetector detector( max_points );
      ip = detect_interest_points(image, detector, max_points);
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    } else if (detector_is_opencv) {

      OpenCvIpDetectorType ocv_type = OPENCV_IP_DETECTOR_TYPE_SIFT;
      if ( interest_operator == "brisk") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_BRISK;
      } else if ( interest_operator == "orb") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_ORB;
      } else if ( interest_operator == "sift") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_SIFT;
      } else if ( interest_operator == "surf") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_SURF;
      }
      OpenCvInterestPointDetector detector(ocv_type, opencv_normalize, describeInDetect, max_points);
      if (has_nodata)
        ip = detect_interest_points(masked_image, detector, max_points);
      else
        ip = detect_interest_points(image, detector, max_points);
    }
#else // End OpenCV section
    } else {
      vw_out() << "Error: Cannot use ORB, BRISK, SIFT, or SURF if not compiled with OpenCV!" << std::endl;
      exit(0); 
    }
#endif

    vw_out() << "Detected " << ip.size() << " raw keypoints!\n";
    if (ip.size() == 0) {
      vw_out() << "No IP, quitting this image!\n";
      continue;
    }

    // Removing Interest Points on nodata or within 1/px
    if (has_nodata) {
      size_t before_size = ip.size();
      remove_ip_near_nodata(raw_image, nodata, ip, nodata_radius );
      
      vw_out(InfoMessage,"interest_point") << "Removed " << before_size-ip.size() << " points close to nodata.\n";
    } // End clearing IP's near nodata

    vw_out() << "\t Found " << ip.size() << " points.\n";

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
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    // Nothing happens here because we compute the descriptors at the same time we detect IPs.
    }
#else
    } else {
      vw_out() << "Error: Cannot use SIFT, SURF, ORB or BRISK if not compiled with OpenCV!" << std::endl;
      exit(0); 
    }
#endif

    // Print raw IP info if requested
    int limit = print_num_ip;
    for (InterestPointList::const_iterator iter=ip.begin(); iter!=ip.end(); ++iter) {
      if (limit <= 0)
        break;
      std::cout << iter->to_string() << "\n";  
      --limit;
    }

    // If ASCII output was requested, write it out.  Otherwise stick with binary output.
    if (vm.count("lowe")) {
      vw_out() << "Writing output file " << file_prefix + ".key" << std::endl;
      write_lowe_ascii_ip_file(file_prefix + ".key", ip);
    } else {
      vw_out() << "Writing output file " << file_prefix + ".vwip" << std::endl;
      write_binary_ip_file(file_prefix + ".vwip", ip);
    }

    // Write Debug image
    if (debug_image > 0) {
      const int WRITE_FULL_RES_DEBUG = 2;
      write_debug_image(file_prefix + "_debug.png",
                        raw_image, ip, has_nodata, nodata,
                        (debug_image==WRITE_FULL_RES_DEBUG),
                        (interest_operator == "orb")); // Orb has large scale values so shrink them.
    }
                         
  } // End loop through input files

  return 0;
}
