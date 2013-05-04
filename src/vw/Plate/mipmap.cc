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


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/PlateManager.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Mipmap takes the following options
// uint32 level
// BBox2i box
// TransactionOrNeg input_transaction
// bool preblur
// ProgressCallback

class MipmapParameters {

  void error(std::string arg, std::string const& params) {
    vw_throw( ArgumentErr() << "Error parsing arguments for --" << arg << " : " << params );
  }

public:
  int starting_level, stopping_level;
  TransactionOrNeg transaction_id;
  BBox2i region;
  std::string mode;

  MipmapParameters( std::string const& mode,
                    std::string const& region_string,
                    std::string const& level_string,
                    TransactionOrNeg tid ) :
    transaction_id(tid), mode(mode) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:@");

    if ( level_string.empty() ) {
      this->error("level", level_string);
    } else {
      tokenizer tokens(level_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("level", level_string);
      starting_level = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("level", level_string);
      stopping_level = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("level", level_string);
    }

    if ( region_string.empty() ) {
      region = BBox2i(0, 0, 0x1<<starting_level, 0x1<<starting_level);
    } else {
      tokenizer tokens(region_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("region", region_string);
      region.min()[0] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      region.min()[1] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      region.max()[0] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      region.max()[1] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("region", region_string);
    }
  }
};

template <class PixelT>
void do_mipmap(boost::shared_ptr<PlateFile> input_plate, MipmapParameters const& mipmap_params) {

  typedef PlateManager<PixelT> PM;
  boost::scoped_ptr<PM> pm(PM::make(mipmap_params.mode, input_plate));

  input_plate->transaction_resume( mipmap_params.transaction_id.promote() );
  input_plate->write_request();

  pm->mipmap( mipmap_params.starting_level, mipmap_params.region,
              mipmap_params.transaction_id, true,
              TerminalProgressCallback("plate.tools.mipmap","mipmap:"),
              mipmap_params.stopping_level);

  input_plate->write_complete();
  input_plate->transaction_end(true);
}

int main( int argc, char *argv[] ) {
  Url plate_url;
  std::string mode, region_string, level_string;
  TransactionOrNeg transaction_id;
  bool help;

  po::options_description general_options("\nUtility for mipmapping a transaction ID that exists only on one level");
  general_options.add_options()
    ("level-range,l", po::value(&level_string), "where arg = <start level>:<stop level>. Start level is where the data currently resides and stop level is where you want the data to finally be at.")
    ("region", po::value(&region_string), "where arg = <ul_x>,<ul_y>:<lr_x>,<lr_y>. Optional.")
    ("transaction-id,t", po::value(&transaction_id), "transaction ID to request and to write.")
    ("mode,m", po::value(&mode), "Output mode [toast, equi, polar]")
    ("help,h", po::bool_switch(&help), "Display this help message.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("url", po::value(&plate_url));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("url", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <url>..." << std::endl << std::endl;
  usage << general_options << std::endl;

  try {
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e ) {
    VW_OUT() << usage.str() << std::endl
             << "Failed to parse command line arguments:" << std::endl
             << "\t" << e.what() << std::endl;
    return 1;
  }

  if ( help || plate_url == Url() ) {
    VW_OUT() << usage.str() << "\n";;
    return 1;
  }

  try {
    MipmapParameters mipmap_params(mode, region_string, level_string, transaction_id);

    VW_ASSERT(mode == "toast" || mode == "equi" || mode == "polar",
              IOErr() << "Unknown projection mode: " << mode);

    boost::shared_ptr<PlateFile> platefile( new PlateFile(plate_url) );
    switch(platefile->pixel_format()) {
    case VW_PIXEL_GRAYA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_mipmap<PixelGrayA<uint8> >(platefile, mipmap_params);
        break;
      case VW_CHANNEL_INT16:
        do_mipmap<PixelGrayA<int16> >(platefile, mipmap_params);
        break;
      case VW_CHANNEL_FLOAT32:
        do_mipmap<PixelGrayA<float32> >(platefile, mipmap_params);
        break;
      default:
        vw_throw(ArgumentErr() << "Image contains a channel type not supported by snapshot.\n");
      }
      break;

    case VW_PIXEL_RGBA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_mipmap<PixelRGBA<uint8> >(platefile, mipmap_params);
        break;
      default:
        vw_throw(ArgumentErr() << "Image contains a channel type not supported by snapshot.\n");
      }
      break;
    default:
      vw_throw(ArgumentErr() << "Image contains a pixel type not supported by snapshot.\n");
    }
  } catch ( const vw::Exception& e ) {
    VW_OUT() << "An error occured: " << e.what() << "\nExiting.\n\n";
    return 1;
  }

  return 0;
}
