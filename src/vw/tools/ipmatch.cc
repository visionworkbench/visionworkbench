// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ipmatch.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///
#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Math.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Camera/CameraGeometry.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

// Draw the two images side by side with matching interest points
// shown with lines.
static void write_match_image(std::string const& out_file_name,
                              std::string const& file1,
                              std::string const& file2,
                              std::vector<InterestPoint> const& matched_ip1,
                              std::vector<InterestPoint> const& matched_ip2) {
  // Skip image pairs with no matches.
  if (matched_ip1.empty())
    return;

  DiskImageView<PixelRGB<uint8> > src1(file1);
  DiskImageView<PixelRGB<uint8> > src2(file2);

  // Work out the scaling to produce the subsampled images. These
  // values are choosen just allow a reasonable rendering time.
  float sub_scale =
    sqrt(1500.0 * 1500.0 / float(src1.cols() * src1.rows()));
  sub_scale +=
    sqrt(1500.0 * 1500.0 / float(src2.cols() * src2.rows()));
  sub_scale /= 2;
  if ( sub_scale > 1 ) sub_scale = 1;

  mosaic::ImageComposite<PixelRGB<uint8> > composite;
  composite.insert(pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(resample(src1.impl(), sub_scale))),
		   0, 0 );
  composite.insert(pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(resample(src2.impl(), sub_scale))),
		   int32(src1.impl().cols() * sub_scale), 0 );
  composite.set_draft_mode( true );
  composite.prepare();

  // Rasterize the composite so that we can draw on it.
  ImageView<PixelRGB<uint8> > comp = composite;

  // Draw a red line between matching interest points
  for (size_t k = 0; k < matched_ip1.size(); ++k) {
    Vector2f start(matched_ip1[k].x, matched_ip1[k].y);
    Vector2f end(matched_ip2[k].x+src1.impl().cols(), matched_ip2[k].y);
    start *= sub_scale;
    end   *= sub_scale;
    float inc_amt = 1/norm_2(end-start);
    for (float r=0; r<1.0; r+=inc_amt ){
      int i = (int)(0.5 + start.x() + r*(end.x()-start.x()));
      int j = (int)(0.5 + start.y() + r*(end.y()-start.y()));
      if (i >=0 && j >=0 && i < comp.cols() && j < comp.rows())
        comp(i,j) = PixelRGB<uint8>(255, 0, 0);
    }
  }

  boost::scoped_ptr<vw::DiskImageResource> rsrc( DiskImageResource::create(out_file_name, comp.format()) );
  block_write_image( *rsrc, comp,
		     TerminalProgressCallback( "tools.ipmatch", "Writing Debug:" ) );
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  double matcher_threshold;
  std::string ransac_constraint;
  float inlier_threshold;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("matcher-threshold,t", po::value(&matcher_threshold)->default_value(0.6), "Threshold for the interest point matcher.")
    ("non-kdtree", "Use an implementation of the interest matcher that is not reliant on a KDTree algorithm")
    ("ransac-constraint,r", po::value(&ransac_constraint)->default_value("similarity"), "RANSAC constraint type.  Choose one of: [similarity, homography, fundamental, or none].")
    ("inlier-threshold,i", po::value(&inlier_threshold)->default_value(10), "RANSAC inlier threshold.")
    ("debug-image,d", "Write out debug images.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value<std::vector<std::string> >(&input_file_names));

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

  if( input_file_names.size() < 2 ) {
    vw_out() << "Error: Must specify at least two input files!" << std::endl << std::endl;
    vw_out() << usage.str();
    return 1;
  }

  // Iterate over combinations of the input files and find interest points in each.
  for (size_t i = 0; i < input_file_names.size(); ++i) {
    for (size_t j = i+1; j < input_file_names.size(); ++j) {

      // Read each file off disk
      std::vector<InterestPoint> ip1, ip2;
      ip1 = read_binary_ip_file(fs::path(input_file_names[i]).replace_extension("vwip").string() );
      ip2 = read_binary_ip_file(fs::path(input_file_names[j]).replace_extension("vwip").string() );
      vw_out() << "Matching between " << input_file_names[i] << " (" << ip1.size() << " points) and " << input_file_names[j] << " (" << ip2.size() << " points).\n";

      std::vector<InterestPoint> matched_ip1, matched_ip2;

      if ( !vm.count("non-kdtree") ) {
        // Run interest point matcher that uses KDTree algorithm.
        InterestPointMatcher<L2NormMetric,NullConstraint> matcher(matcher_threshold);
        matcher(ip1, ip2, matched_ip1, matched_ip2, false,
                TerminalProgressCallback( "tools.ipmatch","Matching:"));
      } else {
        // Run interest point matcher that does not use KDTree algorithm.
        InterestPointMatcherSimple<L2NormMetric,NullConstraint> matcher(matcher_threshold);
        matcher(ip1, ip2, matched_ip1, matched_ip2, false,
                TerminalProgressCallback( "tools.ipmatch","Matching:"));
      }

      remove_duplicates(matched_ip1, matched_ip2);
      vw_out() << "Found " << matched_ip1.size() << " putative matches.\n";

      std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
	ransac_ip2 = iplist_to_vectorlist(matched_ip2);
      std::vector<size_t> indices;
      try {
        // RANSAC is used to fit a transform between the matched sets
        // of points.  Points that don't meet this geometric
        // contstraint are rejected as outliers.
        if (ransac_constraint == "similarity") {
          math::RandomSampleConsensus<math::SimilarityFittingFunctor, math::InterestPointErrorMetric> ransac( math::SimilarityFittingFunctor(),
                                                                                                              math::InterestPointErrorMetric(),
                                                                                                              inlier_threshold ); // inlier_threshold
          Matrix<double> H(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Similarity: " << H << "\n";
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "homography") {
          math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> ransac( math::HomographyFittingFunctor(),
                                                                                                              math::InterestPointErrorMetric(),
                                                                                                              inlier_threshold ); // inlier_threshold
          Matrix<double> H(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Homography: " << H << "\n";
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "fundamental") {
          math::RandomSampleConsensus<camera::FundamentalMatrix8PFittingFunctor, camera::FundamentalMatrixDistanceErrorMetric> ransac( camera::FundamentalMatrix8PFittingFunctor(), camera::FundamentalMatrixDistanceErrorMetric(), inlier_threshold );
          Matrix<double> F(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Fundamental: " << F << "\n";
          indices = ransac.inlier_indices(F,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "none") {
          indices.reserve( matched_ip1.size() );
          for ( size_t i = 0; i < matched_ip1.size(); ++i )
            indices.push_back(i);
        } else {
          vw_out() << "Unknown RANSAC constraint type: " << ransac_constraint
		   << ".  Choose one of: [similarity, homography, fundamental, or none]\n";
	  return 1;
        }
      } catch (const vw::math::RANSACErr& e ) {
        vw_out() << "RANSAC Failed: " << e.what() << "\n";
        continue;
      }
      vw_out() << "Found " << indices.size() << " final matches.\n";

      std::vector<InterestPoint> final_ip1, final_ip2;
      BOOST_FOREACH( size_t& index, indices ) {
        final_ip1.push_back(matched_ip1[index]);
        final_ip2.push_back(matched_ip2[index]);
      }

      std::string output_prefix =
	fs::path(input_file_names[i]).replace_extension().string() + "__" +
        fs::path(input_file_names[j]).stem();
      write_binary_match_file(output_prefix+".match", final_ip1, final_ip2);

      if (vm.count("debug-image")) {
        write_match_image(output_prefix+".tif",
                          input_file_names[i], input_file_names[j],
                          final_ip1, final_ip2);
      }
    }
  }

  return 0;
}

