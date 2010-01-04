// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// Options:
///
///  <none> - If no options are supplied, this utility will create a
///          snapshot of the entire platefile starting with
///          read_cursor and going to the most recent transaction_id.
///
///  --start - Starts a new multi-part snapshot.  Returns the
///          transaction_id for this snapshot as a return value.
///
///  --finish <transaction_id> - Finish a multi-part snapshot. 
///
///   --snapshot <level> <t_begin> <t_end> - Creates a snapshot of the
///          mosaic by compositing together tiles from the
///          transaction_id's in the range [t_begin, t_end] and the
///          level specified (the range is inclusive).
///
///  --region <ul_x> <ul_y> <lr_x> <lr_y> - Limit the snapshot to the
///          region bounded by these upper left (ul) and lower right
///          (lr) coordinates.
///

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ToastPlateManager.h>
//#include <vw/Plate/KmlPlateManager.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::mosaic;
using namespace vw::cartography;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

class SnapshotParameters {
  int m_level;
  int m_min_transaction_id;
  int m_max_transaction_id;
  BBox2i m_region;
  bool m_valid;

  void error(std::string arg, std::string const& params) {
    vw_out(0) << "Error parsing arguments for --" << arg << " : " << params << "\n";
    exit(1);
  }

public:
  SnapshotParameters(std::string const& snapshot_string, std::string const& region_string) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:@");

    if (snapshot_string.size() == 0 && region_string.size() != 0) {
      std::cout << "Error: You must use the --region option in conjunction " 
                << "with the --snapshot option.";
      exit(1);
    }
    
    // STEP 1 : PARSE SNAPSHOT STRING

    if (snapshot_string.size() == 0) {

      // If the argument string is empty, then the user did not specify
      // any snapshot parameters.  We mark ourselves as invalid.
      m_valid = false;

    } else {

      // If the argument string is not empty, we attempt to parse the
      // three parameters out from the snapshot string.
      tokenizer tokens(snapshot_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("snapshot", snapshot_string);
      m_min_transaction_id = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("snapshot", snapshot_string);
      m_max_transaction_id = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("snapshot", snapshot_string);
      m_level = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;
      
      if (tok_iter != tokens.end()) this->error("snapshot", snapshot_string);

      m_valid = true;
    }

    // STEP 2 : PARSE REGION STRING
    
    if (region_string.size() == 0) {

      // If the region string is empty, then the user has not
      // specified the region.  We set the bounding box to the entire
      // level.
      int tiles_per_side = pow(2,m_level);
      m_region = BBox2i(0,0,tiles_per_side,tiles_per_side);

    } else {

      // If the argument string is not empty, we attempt to parse the
      // three parameters out from the snapshot string.
      tokenizer tokens(region_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("region", region_string);
      m_region.min()[0] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      m_region.min()[1] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      m_region.max()[0] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("region", region_string);
      m_region.max()[1] = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("region", region_string);
    }

  }

  bool valid() const { return m_valid; }
  int level() const { return m_level; }
  Vector2i transaction_range() const { return Vector2i(m_min_transaction_id, 
                                                       m_max_transaction_id); }
  BBox2i region() const { return m_region; }
};

// --------------------------------------------------------------------------
//                                DO_MOSAIC
// --------------------------------------------------------------------------

template <class PixelT>
void do_snapshot(boost::shared_ptr<PlateFile> platefile,
                 SnapshotParameters snapshot_parameters) {

  if (snapshot_parameters.valid()) {

    // XXX : Implement me!!

  } else {

    // If the snapshot parameters are not valid, then the user must
    // have run this program without supplying any explicit
    // instructions.  We build a snapshot of the entire platefile
    // starting with read_cursor and going to the most recent
    // transaction_id.
    
    // XXX : Implement me!!

  }
  
  // if (output_mode == "toast") {

  //   boost::shared_ptr<ToastPlateManager<typename ViewT::pixel_type> > pm( 
  //     new ToastPlateManager<typename ViewT::pixel_type> (platefile, num_threads) );

  //   pm->insert(view.impl(), filename, georef, g_debug,
  //              TerminalProgressCallback(InfoMessage, status_str.str()) );

  // }  else if (output_mode == "kml")  {

  //   // boost::shared_ptr<KmlPlateManager> pm = 
  //   //   boost::shared_ptr<KmlPlateManager>( new KmlPlateManager(platefile, num_threads) );

  //   // pm->insert(view.impl(), filename, georef,
  //   //            TerminalProgressCallback(InfoMessage, status_str.str()) );

  // }

}

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  std::string url;
  std::string start_description;
  int finish_transaction_id;
  std::string snapshot_string;
  std::string region_string;

  po::options_description general_options("\nCreate a snapshot of a quadtree.  If no options are supplied, this utility will create a snapshot of the entire platefile starting with read_cursor and going to the max transaction_id.\n\nGeneral options");
  general_options.add_options()
    ("start", po::value<std::string>(&start_description), "where arg = <description> - Starts a new multi-part snapshot.  Returns the transaction_id for this snapshot as a return value.")
    ("finish", po::value<int>(&finish_transaction_id), "where arg = <transaction_id> - Finish a multi-part snapshot.")
    ("snapshot", po::value<std::string>(&snapshot_string), "where arg = <t_begin>:<t_end>@<level> - Creates a snapshot of the mosaic by compositing together tiles from the transaction_id's in the range [t_begin, t_end] at the level specified.") 
    ("region", po::value<std::string>(&region_string), "where arg = <ul_x>,<ul_y>:<lr_x>,<lr_y> - Limit the snapshot to the region bounded by these upper left (ul) and lower right (lr) coordinates.");
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::string>(&url), "");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate_filename> [options]\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try { 
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 1;    
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if (url.size() == 0) {
    std::cout << usage.str();
    return 1;
  }

  // Parse out the rest of the command line options into a special
  // SnapshotParameters class.
  SnapshotParameters snapshot_params(snapshot_string, region_string);

  try {

    //--------------------------- OPEN THE PLATE FILE -----------------------------

    std::cout << "\nOpening plate file: " << url << "\n";
    boost::shared_ptr<PlateFile> platefile = 
      boost::shared_ptr<PlateFile>( new PlateFile(url) );

    //------------------------ START/FINISH TRANSACTION ---------------------------

    if (vm.count("start")) {
      std::vector<TileHeader> dummy;
      int transaction_id = platefile->transaction_request(start_description, dummy);
      vw_out(0) << "Transaction started with ID = " << transaction_id;
      exit(transaction_id);
    }

    if (vm.count("finish")) {
      platefile->transaction_complete(finish_transaction_id);
      vw_out(0) << "Transaction " << finish_transaction_id << " complete.\n";
      exit(0);
    }

    //---------------------------- SNAPSHOT REQUEST -------------------------------
    
    PixelFormatEnum pixel_format = platefile->pixel_format();
    ChannelTypeEnum channel_type = platefile->channel_type();
    
    // Dispatch to the compositer based on the pixel type of this mosaic.
    switch(platefile->pixel_format()) {
    case VW_PIXEL_GRAYA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:  
        do_snapshot<PixelGray<uint8> >(platefile, snapshot_params);
        break;
      default:
        vw_throw(ArgumentErr() << "Image contains a channel type not supported by image2plate.\n");
      }
      break;

    case VW_PIXEL_RGBA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:  
        do_snapshot<PixelGray<uint8> >(platefile, snapshot_params);
        break;
      default:
        std::cout << "Platefile contains a channel type not supported by image2plate.\n";
        exit(0);
      }
      break;
    default:
      std::cout << "Image contains a pixel type not supported by image2plate.\n";
    }

  }  catch (vw::Exception &e) {
    std::cout << "An error occured: " << e.what() << "\nExiting.\n\n";
    exit(1);
  }
     
}
