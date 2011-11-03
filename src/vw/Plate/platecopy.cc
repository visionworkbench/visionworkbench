// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PlateFile.h>
#include <vw/Plate/detail/MipmapHelpers.h>
#include <vw/Plate/TileManipulation.h>
using namespace vw;
using namespace vw::platefile;
namespace d = vw::platefile::detail;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

class CopyParameters {

  void error(std::string arg, std::string const& params) {
    vw_out(ErrorMessage) << "Error parsing arguments for --" << arg << " : " << params << "\n";
    exit(1);
  }

public:
  int level;
  TransactionOrNeg transaction_input_id, transaction_output_id;
  BBox2i region;

  // Constructor
  CopyParameters( std::string const& region_string,
                  TransactionOrNeg input_id,
                  TransactionOrNeg output_id ) :
    transaction_input_id(input_id), transaction_output_id(output_id) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:@");

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

      if (tok_iter == tokens.end()) this->error("region", region_string);
      level = boost::lexical_cast<int>(*tok_iter);
      ++tok_iter;

      if (tok_iter != tokens.end()) this->error("region", region_string);
    }
  }
};

template <class PixelT>
void copy_job( int level, BBox2i const& region,
               TransactionOrNeg input_tid,
               boost::shared_ptr<ReadOnlyPlateFile> input_plate,
               boost::shared_ptr<PlateFile> output_plate,
               const ProgressCallback &progress ) {
  const uint64 TILE_BYTES  = input_plate->default_tile_size() * input_plate->default_tile_size() * uint32(PixelNumBytes<PixelT>::value);
  const uint64 CACHE_BYTES = vw_settings().system_cache_size();
  const uint64 CACHE_TILES = CACHE_BYTES / TILE_BYTES;
  VW_ASSERT( CACHE_TILES > 1, LogicErr() << "Cache too small to process any tiles in snapshot." );
  // This is an arbitrary value, to hopefully catch
  // pathlogically-small cache sizes
  if ( CACHE_TILES < 100 )
    vw_out(WarningMessage) << "You will lose a lot speed to thrashing if you can't cache at least 100 tiles (you can only store " << CACHE_TILES << ")\n";

  std::list<BBox2i> sub_regions = bbox_tiles(region, 512, 512);
  float inc_amt = 1./float(sub_regions.size());
  progress.report_progress(0);
  BOOST_FOREACH( const BBox2i& sub_region, sub_regions ) {
    std::list<TileHeader> hdrs = input_plate->search_by_region( level, sub_region, input_tid );
    size_t hdrs_size = hdrs.size();
    for ( size_t i = 0; i < hdrs_size; i += CACHE_TILES ) {
      std::list<TileHeader>::const_iterator hdrs_start = hdrs.begin(), hdrs_stop = hdrs.begin();
      std::advance(hdrs_start, i);
      std::advance(hdrs_stop, i + ((i+CACHE_TILES < hdrs_size) ? CACHE_TILES : hdrs_size - i ) );
      Datastore::TileSearch tile_lookup;
      std::copy( hdrs_start, hdrs_stop, std::back_inserter(tile_lookup));
      // Load tiles from input and immediately write to output with
      // out decoding the imagery.
      BOOST_FOREACH( const Tile& t, input_plate->batch_read(tile_lookup) ) {
        output_plate->write_update( &t.data->operator[](0), t.data->size(), t.hdr.col(), t.hdr.row(), level, t.hdr.filetype() );
      }
    }
    progress.report_incremental_progress(inc_amt);
  }
  progress.report_finished();
}

template <class PixelT>
void do_copy(boost::shared_ptr<ReadOnlyPlateFile> input_plate,
             boost::shared_ptr<PlateFile> output_plate,
             CopyParameters& copy_parameters ) {
  if (copy_parameters.level != -1 ) {
    output_plate->transaction_resume(copy_parameters.transaction_output_id.promote());
    output_plate->audit_log()
      << "Started multi-part copy (t_id = " << copy_parameters.transaction_output_id
      << ") -- level:" << copy_parameters.level
      << " region: " << copy_parameters.region << "\n";
    output_plate->write_request();

    copy_job<PixelT>(copy_parameters.level, copy_parameters.region,
                     copy_parameters.transaction_input_id,
                     input_plate, output_plate,
                     TerminalProgressCallback("plate.tools.platecopy",""));

    output_plate->write_complete();
    output_plate->audit_log()
      << "Finish multi-part copy (t_id = " << copy_parameters.transaction_output_id
      << ") -- level:" << copy_parameters.level
      << " region: " << copy_parameters.region << "\n";
  } else {
    // If no region was specified, then copy the entire transaction
    output_plate->transaction_begin("Full copy (requested t_if: " + vw::stringify(copy_parameters.transaction_output_id) + ")",
                                    copy_parameters.transaction_output_id );

    for ( int32 level = 0; level < input_plate->num_levels(); level++ ) {
      copy_job<PixelT>(level, d::move_down(BBox2i(0,0,1,1), level),
                       copy_parameters.transaction_input_id,
                       input_plate, output_plate,
                       TerminalProgressCallback("plate.tools.platecopy","Level: " + stringify(level)));
    }

    output_plate->transaction_end(true);
  }
}

int main( int argc, char *argv[] ) {

  Url input_url, output_url;
  std::string start_description;
  std::string region_string;
  TransactionOrNeg transaction_input_id = -1, transaction_output_id = -1;

  po::options_description general_options("Copies transactions from one plate to another.");
  general_options.add_options()
    ("start", po::value(&start_description),  "where arg = <description> - Starts a new multi-part snapshot.  Returns the transaction_id for this snapshot as a return value.")
    ("finish", "Finish a multi-part snapshot.")
    ("transaction-output,t",  po::value(&transaction_output_id), "Transaction ID to write to output plate")
    ("transaction-input,i", po::value(&transaction_input_id), "Transaction ID to read from input plate")
    ("region", po::value(&region_string), "where arg = <ul_x>, <ul_y>:<lr_x>, <lr_y>@<level> - Limit the snapshot to the region bounded by these upper left (ul) and lower right (lr) coordinates at the level specified.")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-url", po::value(&input_url), "")
    ("output-url", po::value(&output_url), "");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-url", 1);
  p.add("output-url", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <input plate> <output plate> [options]\n";
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

  if( vm.count("help") || input_url == Url() || output_url == Url() ) {
    std::cout << usage.str();
    return 1;
  }

  try {
    boost::shared_ptr<ReadOnlyPlateFile> input_plate;
    boost::shared_ptr<PlateFile> output_plate;

    VW_ASSERT( output_url != input_url,
               IOErr() << "Input and output url must be different URLs." );
    input_plate.reset( new ReadOnlyPlateFile( input_url ) );
    output_plate.reset( new PlateFile( output_url, "", "", input_plate->default_tile_size(),
                                       input_plate->default_file_type(),
                                       input_plate->pixel_format(), input_plate->channel_type() ) );

    if (vm.count("start")) {
      Transaction t = output_plate->transaction_begin(start_description, transaction_output_id);
      vw_out() << "Transaction started with ID = " << t << "\n";
      exit(0);
    }

    if (vm.count("finish")) {
      output_plate->transaction_resume(transaction_output_id.promote());
      output_plate->transaction_end(true);
      vw_out() << "Transaction " << transaction_output_id << " complete.\n";
      exit(0);
    }

    CopyParameters copy_params(region_string,
                               transaction_input_id,
                               transaction_output_id);

    // Dispatch to the copy machine
    switch(input_plate->pixel_format()) {
    case VW_PIXEL_GRAYA:
      switch(input_plate->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_copy<PixelGrayA<uint8> >(input_plate, output_plate, copy_params); break;
      case VW_CHANNEL_INT16:
        do_copy<PixelGrayA<int16> >(input_plate, output_plate, copy_params); break;
      case VW_CHANNEL_FLOAT32:
        do_copy<PixelGrayA<float32> >(input_plate, output_plate, copy_params); break;
      default:
        vw_throw(ArgumentErr() << "Plate contains a channel type not supported by platecopy.\n");
        exit(1);
      }
      break;
    case VW_PIXEL_RGBA:
      switch(input_plate->channel_type()) {
      case VW_CHANNEL_UINT8:
        do_copy<PixelRGBA<uint8> >(input_plate, output_plate, copy_params); break;
      default:
        vw_throw(ArgumentErr() << "Plate contains a channel type not supported by platecopy.\n");
        exit(1);
      }
      break;
    default:
      std::cout << "Plate contains a pixel type not supported by snapshot.\n";
      exit(1);
    }
  } catch ( const vw::Exception& e ) {
    std::cout << "An error occured: " << e.what() << "\nExiting.\n\n";
    exit(1);
  }

  return 0;
}
