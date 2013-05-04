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
    vw_throw( ArgumentErr() << "Error parsing arguments for --" << arg << " : " << params );
  }

public:

  // Public variables.  You access these directly.
  int level;
  TransactionOrNeg begin_transaction_id;
  TransactionOrNeg end_transaction_id;
  TransactionOrNeg write_transaction_id;
  bool tweak_settings_for_terrain;
  BBox2i region;

  // Constructor
  SnapshotParameters(std::string const& range_string,
                     std::string const& region_string,
                     TransactionOrNeg write_transaction_id,
                     bool terrain) :
    write_transaction_id(write_transaction_id),
    tweak_settings_for_terrain(terrain) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:@");

    // STEP 1 : PARSE RANGE STRING

    if (range_string.empty()) {
      vw_throw( ArgumentErr() << "Error: You must specify a transaction_id range "
                << "with the --transaction-range option.\n" );
    } else {

      // If the range string is not empty, we attempt to parse the
      // two parameters out from the range string.
      tokenizer tokens(range_string, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (tok_iter == tokens.end()) this->error("snapshot", range_string);
      begin_transaction_id = boost::lexical_cast<TransactionOrNeg>(*tok_iter);
      ++tok_iter;

      if (tok_iter == tokens.end()) this->error("snapshot", range_string);
      end_transaction_id = boost::lexical_cast<TransactionOrNeg>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("snapshot", range_string);

    }

    // STEP 2 : PARSE REGION STRING

    if (region_string.empty()) {

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
void do_snapshot(boost::shared_ptr<ReadOnlyPlateFile> input_plate, boost::shared_ptr<PlateFile> output_plate, SnapshotParameters snapshot_parameters) {

  SnapshotManager<PixelT> sm( input_plate, output_plate, snapshot_parameters.tweak_settings_for_terrain );
  if (snapshot_parameters.level != -1) {

    if (snapshot_parameters.write_transaction_id == -1) {
      vw_throw( ArgumentErr() << "Error: you must specify a transaction_id for this snapshot "
                << "using the --transaction_id flag." );
    }

    output_plate->audit_log()
      << "Started multi-part snapshot (t_id = " << snapshot_parameters.write_transaction_id
      << ") -- level:" << snapshot_parameters.level
      << " region: " << snapshot_parameters.region << "\n";

    // Grab a lock on a blob file to use for writing tiles during
    // the two operations below.
    output_plate->write_request();

    // If the user has specified a region, then we
    sm.snapshot(snapshot_parameters.level, snapshot_parameters.region,
                TransactionRange( snapshot_parameters.begin_transaction_id, snapshot_parameters.end_transaction_id));

    // Release the blob id lock.
    output_plate->write_complete();

    output_plate->audit_log()
      << "Finished multi-part snapshot (t_id = " << snapshot_parameters.write_transaction_id
      << ") -- level:" << snapshot_parameters.level
      << " region: " << snapshot_parameters.region << "\n";

  } else {

    // If no region was specified, then we build a snapshot of the entire
    // output_plate.

    // Request a transaction (write_transaction_id might be -1, and that's okay)
    const std::string description = "Full snapshot (requested t_id:  " + vw::stringify(snapshot_parameters.write_transaction_id) + ")";
    output_plate->transaction_begin(description, snapshot_parameters.write_transaction_id);

    // Grab a lock on a blob file to use for writing tiles during
    // the two operations below.
    output_plate->write_request();

    // Do a full snapshot
    sm.full_snapshot(TransactionRange(snapshot_parameters.begin_transaction_id, snapshot_parameters.end_transaction_id));

    // Release the blob id lock and note that the transaction is finished.
    output_plate->write_complete();
    output_plate->transaction_end(true);
  }

}

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  Url input_url, output_url;
  std::string start_description;
  std::string output_mode;
  std::string range_string;
  std::string region_string;
  TransactionOrNeg transaction_id = -1;

  po::options_description general_options("\nCreate a snapshot of a quadtree.  If no options are supplied, this utility will create a snapshot of the entire platefile starting with read_cursor and going to the max transaction_id.\n\nGeneral options");
  general_options.add_options()
    ("start", po::value(&start_description) , "where arg = <description> - Starts a new multi-part snapshot.  Returns the transaction_id for this snapshot as a return value.")
    ("finish", "Finish a multi-part snapshot.")
    ("transaction-id,t",  po::value(&transaction_id)                                                                                                                     , "Transaction ID to use for starting/finishing/or snapshotting.")
    ("transaction-range", po::value(&range_string)      , "where arg = <t_begin>:<t_end> - Creates a snapshot of the mosaic by compositing together tiles from the transaction_ids in the range [t_begin , t_end] (inclusive).")
    ("region",       po::value(&region_string)     , "where arg = <ul_x>                                                                                                                            , <ul_y>:<lr_x>                                                    , <lr_y>@<level> - Limit the snapshot to the region bounded by these upper left (ul) and lower right (lr) coordinates at the level specified.")
    ("terrain",      "Tweak settings for terrain.")
    ("output-plate,o", po::value(&output_url),         "output to this plate instead of the input plate")
    ("help,h",       "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value(&input_url), "");

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
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if (vm.count("input-file") != 1) {
    std::cout << usage.str();
    return 1;
  }

  try {

    //--------------------------- OPEN THE PLATE FILE -----------------------------
    boost::shared_ptr<ReadOnlyPlateFile> input_plate;
    boost::shared_ptr<PlateFile> output_plate;

    if (output_url != Url() && (output_url != input_url)) {
      std::cout << "\nOpening input plate file: " << input_url.string() << "\n";
      input_plate.reset(new ReadOnlyPlateFile(input_url));
      std::cout << "\nOpening output plate file: " << output_url.string() << "\n";
      output_plate.reset(new PlateFile(output_url, input_plate->index_header().type(), input_plate->index_header().description(), input_plate->index_header().tile_size(), input_plate->index_header().tile_filetype(), input_plate->pixel_format(), input_plate->channel_type()));
    } else {
      std::cout << "\nOpening plate file: " << input_url.string() << "\n";
      output_plate.reset(new PlateFile(input_url));
      input_plate = output_plate;
    }

    //------------------------ START/FINISH TRANSACTION ---------------------------

    if (vm.count("start")) {
      if (!vm.count("transaction-id")) {
        VW_OUT() << "You must specify a transaction-id if you use --start.\n";
        return 1;
      }

      Transaction t = output_plate->transaction_begin(start_description, transaction_id);
      VW_OUT() << "Transaction started with ID = " << t << "\n";
      VW_OUT() << "Plate has " << input_plate->num_levels() << " levels.\n";
      return 0;
    }

    if (transaction_id.newest()) {
      VW_OUT() << "You must specify a transaction-id if you are not starting a new one.\n";
      return 1;
    }

    output_plate->transaction_resume(transaction_id.promote());

    if (vm.count("finish")) {
      // Update the read cursor when the snapshot is complete!
      output_plate->transaction_end(true);
      VW_OUT() << "Transaction " << transaction_id << " complete.\n";
      return 0;
    }

    //---------------------------- SNAPSHOT REQUEST -------------------------------

    // Parse out the rest of the command line options into a special
    // SnapshotParameters class.
    SnapshotParameters snapshot_params(range_string, region_string,
                                       transaction_id, vm.count("terrain"));

    // Dispatch to the compositer based on the pixel type of this mosaic.
    switch(input_plate->pixel_format()) {
    case VW_PIXEL_GRAYA:
      switch(input_plate->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_snapshot<PixelGrayA<uint8> >(input_plate, output_plate, snapshot_params);
        break;
      case VW_CHANNEL_INT16:
        do_snapshot<PixelGrayA<int16> >(input_plate, output_plate, snapshot_params);
        break;
      case VW_CHANNEL_FLOAT32:
        do_snapshot<PixelGrayA<float32> >(input_plate, output_plate, snapshot_params);
        break;
      default:
        vw_throw(ArgumentErr() << "Image contains a channel type not supported by snapshot.\n");
        return 1;
      }
      break;

    case VW_PIXEL_RGBA:
      switch(input_plate->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_snapshot<PixelRGBA<uint8> >(input_plate, output_plate, snapshot_params);
        break;
      default:
        std::cout << "Platefile contains a channel type not supported by snapshot.\n";
        return 1;
      }
      break;
    default:
      std::cout << "Image contains a pixel type not supported by snapshot.\n";
      return 1;
    }

  }  catch (const vw::Exception& e) {
    std::cout << "An error occured: " << e.what() << "\nExiting.\n\n";
    return 1;
  }

  return 0;
}
