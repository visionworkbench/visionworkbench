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
///          In this mode, the snapshot tool starts and finishes a
///          transaction ID automatically.
///
///  --start - Starts a new multi-part snapshot.  Returns the
///          transaction_id for this snapshot as a return value.
///
///  --finish <transaction_id> - Finish a multi-part snapshot. 
///
///   --snapshot <t_begin>:<t_end>:<t_write>@<level> - Creates a
///          snapshot of the mosaic by compositing together tiles from
///          the transaction_id's in the range [t_begin, t_end] and
///          the level specified (the range is inclusive).  New tiles
///          in the snapshot are written using the t_write
///          transaction_id.
///
///  --region <ul_x>,<ul_y>:<lr_x>,<lr_y> - Limit the snapshot to the
///          region bounded by these upper left (ul) and lower right
///          (lr) coordinates.
///

#include <vw/Plate/PlateFile.h>
#include <vw/Plate/SnapshotManager.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

class SnapshotParameters {

  void error(std::string arg, std::string const& params) {
    vw_out(ErrorMessage) << "Error parsing arguments for --" << arg << " : " << params << "\n";
    exit(1);
  }

public:

  // Public variables.  You access these directly.
  int level;
  int begin_transaction_id;
  int end_transaction_id;
  int write_transaction_id;
  BBox2i region;

  // Constructor
  SnapshotParameters(std::string const& range_string, 
                     std::string const& region_string, 
                     int const& write_transaction_id) : write_transaction_id(write_transaction_id) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:@");

    // STEP 1 : PARSE RANGE STRING

    if (range_string.size() == 0) {
      std::cout << "Error: You must specify a transaction_id range "
                << "with the --range option.\n";
      exit(1);

    } else {

      // If the range string is not empty, we attempt to parse the
      // two parameters out from the range string.
      tokenizer tokens(range_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("snapshot", range_string);
      begin_transaction_id = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("snapshot", range_string);
      end_transaction_id = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;
      
      if (tok_iter != tokens.end()) this->error("snapshot", range_string);

    }

    // STEP 2 : PARSE REGION STRING
    
    if (region_string.size() == 0) {

      // If the region string is empty, then the user has not
      // specified the region.  We record this by setting the level to
      // -1;
      level = -1;

    } else {

      // If the argument string is not empty, we attempt to parse the
      // three parameters out from the snapshot string.
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

      if (tok_iter == tokens.end()) this->error("snapshot", region_string);
      level = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("region", region_string);
    }

  }
};

// --------------------------------------------------------------------------
//                                DO_MOSAIC
// --------------------------------------------------------------------------

template <class PixelT>
void do_snapshot(boost::shared_ptr<PlateFile> platefile,
                 SnapshotParameters snapshot_parameters) {

  SnapshotManager<PixelT> sm( platefile );
  if (snapshot_parameters.level != -1) {

    if (snapshot_parameters.write_transaction_id == -1) {
      vw_out(ErrorMessage) << "Error: you must specify a transaction_id for this snapshot "
                           << "using the --transaction_id flag.";
      exit(1);
    }
    
    // Grab a lock on a blob file to use for writing tiles during
    // the two operations below.
    platefile->write_request(0);
      
    // If the user has specified a region, then we 
    sm.snapshot(snapshot_parameters.level, snapshot_parameters.region,
                snapshot_parameters.begin_transaction_id,
                snapshot_parameters.end_transaction_id,
                snapshot_parameters.write_transaction_id);
    
    // Release the blob id lock.
    platefile->write_complete();

  } else {

    // If no region was specified, then we build a snapshot of the
    // entire platefile.
    if (snapshot_parameters.write_transaction_id == -1) {

      // User did not supply a t_id.  We must request and complete a
      // transaction on our own.
      int t_id = platefile->transaction_request("Full snapshot.", -1);

      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      platefile->write_request(0);

      // Do a full snapshot
      sm.full_snapshot(snapshot_parameters.begin_transaction_id,
                       snapshot_parameters.end_transaction_id,
                       t_id);

      // Release the blob id lock and note that the transaction is finished.
      platefile->write_complete();
      platefile->transaction_complete(t_id, true);

    } else {
      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      platefile->write_request(0);

      // User-supplied transaction_id
      sm.full_snapshot(snapshot_parameters.begin_transaction_id,
                       snapshot_parameters.end_transaction_id,
                       snapshot_parameters.write_transaction_id);

      // Release the blob id lock and note that the transaction is finished.
      platefile->write_complete();
    }
  }
           
}

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  std::string url;
  std::string start_description;
  std::string output_mode;
  std::string range_string;
  std::string region_string;
  int transaction_id = -1;

  po::options_description general_options("\nCreate a snapshot of a quadtree.  If no options are supplied, this utility will create a snapshot of the entire platefile starting with read_cursor and going to the max transaction_id.\n\nGeneral options");
  general_options.add_options()
    ("start", po::value<std::string>(&start_description), "where arg = <description> - Starts a new multi-part snapshot.  Returns the transaction_id for this snapshot as a return value.")
    ("finish", "Finish a multi-part snapshot.")
    ("transaction-id,t", po::value<int>(&transaction_id), "Transaction ID to use for starting/finishing/or snapshotting.")
    ("transaction-range", po::value<std::string>(&range_string), "where arg = <t_begin>:<t_end> - Creates a snapshot of the mosaic by compositing together tiles from the transaction_ids in the range [t_begin, t_end] (inclusive).") 
    ("region", po::value<std::string>(&region_string), "where arg = <ul_x>,<ul_y>:<lr_x>,<lr_y>@<level> - Limit the snapshot to the region bounded by these upper left (ul) and lower right (lr) coordinates at the level specified.")
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

  try {

    //--------------------------- OPEN THE PLATE FILE -----------------------------

    std::cout << "\nOpening plate file: " << url << "\n";
    boost::shared_ptr<PlateFile> platefile = 
      boost::shared_ptr<PlateFile>( new PlateFile(url) );

    //------------------------ START/FINISH TRANSACTION ---------------------------

    if (vm.count("start")) {
      if (!vm.count("transaction-id")) {
        std::cout << "You must specify a transaction-id if you use --start.\n";
        exit(1);
      }
      int t = platefile->transaction_request(start_description, transaction_id);
      vw_out() << "Transaction started with ID = " << t << "\n";
      vw_out() << "Plate has " << platefile->num_levels() << " levels.\n";
      exit(0);
    }

    if (vm.count("finish")) {
      if (!vm.count("transaction-id")) {
        vw_out() << "You must specify a transaction-id if you use --finish.\n";
        exit(1);
      }
      // Update the read cursor when the snapshot is complete!
      platefile->transaction_complete(transaction_id, true);
      vw_out() << "Transaction " << transaction_id << " complete.\n";
      exit(0);
    }

    //---------------------------- SNAPSHOT REQUEST -------------------------------
    
    // Parse out the rest of the command line options into a special
    // SnapshotParameters class.
    SnapshotParameters snapshot_params(range_string, region_string, transaction_id);
    PixelFormatEnum pixel_format = platefile->pixel_format();
    ChannelTypeEnum channel_type = platefile->channel_type();
    
    // Dispatch to the compositer based on the pixel type of this mosaic.
    switch(platefile->pixel_format()) {
    case VW_PIXEL_GRAYA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:  
        do_snapshot<PixelGrayA<uint8> >(platefile, snapshot_params);
        break;
      default:
        vw_throw(ArgumentErr() << "Image contains a channel type not supported by image2plate.\n");
      }
      break;

    case VW_PIXEL_RGBA:
      switch(platefile->channel_type()) {
      case VW_CHANNEL_UINT8:  
        do_snapshot<PixelRGBA<uint8> >(platefile, snapshot_params);
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
