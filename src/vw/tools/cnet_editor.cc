// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// boost
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// standard
#include <vector>
#include <iostream>

// VisionWorkbench
#include <vw/Math.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/InterestPoint/InterestData.h>

using namespace vw;
using namespace vw::ip;
using namespace vw::ba;

#if VW_BOOST_VERSION < 103400
std::istream& operator>>(std::istream& is, fs::path& p) {
    std::string s;
    is >> s;
    p = s;
    return is;
}
#endif

// Main Executable
int main( int argc, char *argv[] ) {
  std::string cnet_file;
  std::string image_mean_file;
  fs::path output_cnet_file;
  fs::path output_deleted_index_file;
  ControlNetwork cnet("ControlNetwork Editor");
  double cutoff_sigma;
  std::string output_cnet_file_tmp;
  std::string output_deleted_index_file_tmp;
  std::string data_dir;

  po::options_description general_options("Options");
  general_options.add_options()
    ("cutoff_sigma,c", po::value<double>(&cutoff_sigma)->default_value(2),"This is the highend cutoff sigma.")
    ("output_cnet_file,o", po::value<std::string>(&output_cnet_file_tmp)->default_value("processed.cnet"), "Name of processed control network file to write out.")
    ("output_deleted_index_file,o", po::value<std::string>(&output_deleted_index_file_tmp)->default_value("cnet_deleted.txt"), "Name of file of indexes to deleted control points to write out.")
    ("data_dir,d",po::value<std::string>(&data_dir)->default_value("."), "Name of directory to write output files into.")
    ("clip-if-not-on-moon", "Clip control points that don't track features that are not at least 200 km near the surface of the moon.")
    ("help,h","Brings up this.");

  po::options_description positional_options("Positional Options");
  positional_options.add_options()
    ("cnet", po::value<std::string>(&cnet_file), "Control Network")
    ("image-mean", po::value<std::string>(&image_mean_file), "Image Mean file");

  po::positional_options_description positional_options_desc;
  positional_options_desc.add("cnet", 1);
  positional_options_desc.add("image-mean",1);

  po::options_description all_options("Allowed Options");
  all_options.add(general_options).add(positional_options);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <cnet> <image-mean>\n\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_options_desc).run(), vm );
    po::notify( vm );
  } catch (po::error &e ) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if ( vm.count("help") ||
       !vm.count("cnet") ||
       !vm.count("image-mean") ) {
    vw_out() << usage.str() << std::endl;
    return 1;
  }

  output_cnet_file = output_cnet_file_tmp;
  output_deleted_index_file = output_deleted_index_file_tmp;

  // Loading control network file
  std::vector<std::string> tokens;
  boost::split( tokens, cnet_file, boost::is_any_of(".") );
  if ( tokens[tokens.size()-1] == "net" ) {
    cnet.read_isis( cnet_file );
  } else if ( tokens[tokens.size()-1] == "cnet" ) {
    cnet.read_binary( cnet_file );
  } else {
    vw_throw( IOErr() << "Unknown Control Network file extension, \""
                << tokens[tokens.size()-1] << "\"." );
  }

  // Loading image mean file
  std::ifstream f;
  f.open( image_mean_file.c_str(), std::ios::binary | std::ios::in );
  unsigned read_cnet_size;
  f.read((char*)&(read_cnet_size), sizeof(unsigned));
  if ( read_cnet_size != cnet.size() )
    vw_throw( IOErr() << "Image mean file doesn't seem to match Cnet" );
  unsigned error_size;
  std::list<double> image_errors;
  f.read((char*)&(error_size), sizeof(unsigned));

  for ( unsigned i = 0; i < error_size; i++ ) {
    double temp;
    f.read((char*)&(temp), sizeof(double));
    image_errors.push_back(temp);
  }

  // Collecting statistics
  double min_image = *(std::min_element(image_errors.begin(),
                                       image_errors.end()));
  double max_image = *(std::max_element(image_errors.begin(),
                                        image_errors.end()));
  double list_size = image_errors.size();
  double mean_image=0, stddev_image=0;
  for( std::list<double>::iterator it = image_errors.begin();
       it != image_errors.end(); it++ ) {
    mean_image += (*it) / error_size;
    stddev_image += (*it)*(*it) / error_size;
  }
  stddev_image = sqrt( stddev_image - mean_image*mean_image );
  vw_out() << "Image min: " << min_image << " max: " << max_image
            << " mean: " << mean_image << " stddev: " << stddev_image
            << std::endl;

  // Awesome, now clipping based one std_dev
  int clipping_count = 0;
  int cp_clip_count = 0;
  std::list<double>::iterator image_error = image_errors.begin();
  float inc_amt = 1.0/float(cnet.size());
  TerminalProgressCallback tpc("","Clipping");
  for ( unsigned cpi = 0; cpi < cnet.size(); cpi++) {
    tpc.report_incremental_progress(inc_amt);
    for ( unsigned cmi = 0; cmi < cnet[cpi].size(); cmi++ ) {
      if ( image_error == image_errors.end() )
        vw_throw( IOErr() << "Internal overflow error" );

      // Do clipping
      if ( *image_error >= mean_image + stddev_image*cutoff_sigma ) {
        clipping_count++;

        cnet[cpi].delete_measure( cmi );
        cmi--;
      }
      image_error++;
    }

    if ( cnet[cpi].size() < 2 ||
         ( vm.count("clip-if-not-on-moon") &&
           (fabs(norm_2(cnet[cpi].position())-1737.4e3) > 200e3 ) ) ) {
      clipping_count += int(cnet[cpi].size());
      cnet.delete_control_point( cpi );
      cpi--;
      cp_clip_count++;
    }
  }
  tpc.report_finished();
  if ( image_error != image_errors.end() )
    vw_throw( IOErr() << "Internal overflow error" );
  vw_out() << clipping_count << " control measures removed.\n";
  vw_out() << cp_clip_count << " of control points removed.\n";

  vw_out() << "\nWriting out new control network\n";
  std::string outfile_str = fs::path(data_dir / output_cnet_file ).string();
  vw_out() << "\tfile: " << outfile_str << "\n";
  cnet.write_binary(outfile_str);
}
