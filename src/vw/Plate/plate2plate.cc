// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastDem.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Core/Debugging.h>
#include <boost/foreach.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>

using namespace std;
using namespace vw;
using namespace vw::platefile;

// There isn't much abstraction here. A Filter takes a platefile and writes a
// new platefile. This rules out lazy filters for the moment. The
// col/row/level/transaction_id is for the input value, the output can feel
// free to write things elsewhere.

// The function will be called for every input tile, including all transaction
// ids. The output function should feel no obligation to write a tile for every
// input.

template <typename ImplT>
struct FilterBase {
  inline ImplT& impl()             { return static_cast<ImplT&>(*this); }
  inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

  inline void init(PlateFile& output, const PlateFile& input,
                   TransactionOrNeg input_transaction_id, Transaction output_transaction_id) {
    return impl().init(output, input, input_transaction_id, output_transaction_id);
  }
  inline void fini(PlateFile& output, const PlateFile& input, Transaction output_transaction_id) {
    return impl().fini(output, input, output_transaction_id);
  }

#   define lookup(name, type) type name(type data) const { return impl().name(data); }
    lookup(mode, string);
    lookup(tile_size, int);
    lookup(filetype, string);
    lookup(pixel_format, PixelFormatEnum);
    lookup(channel_type, ChannelTypeEnum);
#   undef lookup

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, TransactionOrNeg input_transaction_id, Transaction output_transaction_id) {
    impl()(output, input, col, row, level, input_transaction_id, output_transaction_id);
  }
};

struct Identity : public FilterBase<Identity> {
#   define lookup(name, type) type name(type data) const { return data; }
    lookup(mode, string);
    lookup(tile_size, int);
    lookup(filetype, string);
    lookup(pixel_format, PixelFormatEnum);
    lookup(channel_type, ChannelTypeEnum);
#   undef lookup

  inline void init(PlateFile& output, const PlateFile& /*input*/, TransactionOrNeg /* input_transaction_id */, Transaction /*output_transaction_id*/) { output.write_request(); }
  inline void fini(PlateFile& output, const PlateFile& /*input*/, Transaction /*transaction_id*/) { output.write_complete(); }

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, TransactionOrNeg input_transaction_id, Transaction output_transaction_id) {
    ImageView<PixelRGBA<double> > tile;
    TileHeader hdr = input.read(tile, col, row, level, input_transaction_id);
    output.write_update(tile, col, row, level, output_transaction_id);
  }
};

struct ToastDem : public FilterBase<ToastDem> {
  string mode(string) const { return "toast_dem"; }
  int tile_size(int)  const { return 32; }
  string filetype(string) const { return "toast_dem_v1"; }
  PixelFormatEnum pixel_format(PixelFormatEnum) const { return VW_PIXEL_SCALAR; }
  ChannelTypeEnum channel_type(ChannelTypeEnum) const { return VW_CHANNEL_UINT8; }

  inline void init(PlateFile& output, const PlateFile& input, TransactionOrNeg input_transaction_id, Transaction output_transaction_id) {
    output.write_request();

    // Write null tiles for the levels we don't have data for
    int level_difference = log(input.default_tile_size()/float(output.default_tile_size())) / log(2.) + 0.5;

    vw_out() << "Creating null tiles for a level difference of " << level_difference << std::endl;

    uint64 bytes;
    boost::shared_array<uint8> null_tile = toast_dem_null_tile(bytes);

    for (int level = 0; level < level_difference; ++level) {
      int region_size = 1 << level;
      for (int row = 0; row < region_size; ++row) {
        for (int col = 0; col < region_size; ++col) {
          DemWriter writer(output);
          make_toast_dem_tile(writer, input, col, row, level, 0,
                              input_transaction_id, output_transaction_id);
        }
      }
    }
  }

  inline void fini(PlateFile& output, const PlateFile& /*input*/, Transaction /*transaction_id*/) {
    output.write_complete();
  }

  struct DemWriter : public ToastDemWriter {
    PlateFile& platefile;
    DemWriter(PlateFile& output) : platefile(output) { }
    inline void operator()(const boost::shared_array<uint8> data, size_t data_size,
                           int32 dem_col, int32 dem_row,
                           int32 dem_level, Transaction output_transaction_id) const {
      platefile.write_update(data, data_size, dem_col, dem_row,
                             dem_level, output_transaction_id);
    }
  };

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, TransactionOrNeg input_transaction_id, Transaction output_transaction_id) {
    DemWriter writer(output);
    int level_difference = log(input.default_tile_size()/
                               float(output.default_tile_size())) / log(2.) + 0.5;
    make_toast_dem_tile(writer, input, col, row, level, level_difference,
                        input_transaction_id, output_transaction_id);
  }
};

struct Options {
  string input_name;
  string output_name;

  string mode;
  string description;
  int tile_size;
  string filetype;
  PixelFormatEnum pixel_format;
  ChannelTypeEnum channel_type;
  int bottom_level;
  bool skim_mode;

  string filter;

  Options() :
    tile_size(0), pixel_format(VW_PIXEL_UNKNOWN), channel_type(VW_CHANNEL_UNKNOWN) {}
};

VW_DEFINE_EXCEPTION(Usage, Exception);

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description options("Options");
  options.add_options()
    ("output-name,o",    po::value(&opt.output_name),  "Specify the URL of the output platefile.")
    ("input-name,i",     po::value(&opt.input_name),   "Specify the URL of the input platefile.")
    ("file-type",        po::value(&opt.filetype),     "Output file type")
    ("mode",             po::value(&opt.mode),         "Output mode [toast, kml]")
    ("tile-size",        po::value(&opt.tile_size),    "Output size, in pixels")
    ("filter",           po::value(&opt.filter),       "Filters to run [identity, toast_dem]")
    ("bottom-level",     po::value(&opt.bottom_level), "Bottom level to process")
    ("skim-last-id-only", "Only process the last transaction id from the input")
    ("help,h",           "Display this help message.");

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw(Usage() << "Error parsing input:\n\t" << e.what() << "\n" << options );
  }

  opt.skim_mode = vm.count("skim-last-id-only");

  if (opt.output_name.empty() || opt.input_name.empty())
    vw_throw(Usage() << "Requires input and output defined!\n" << options);
  if ( opt.filter.empty() )
    vw_throw(Usage() << "Requires filter to be defined!\n" << options);
}

template <typename FilterT>
void run(Options& opt, FilterBase<FilterT>& filter) {
  // XXX: input_name should probably be a url
  PlateFile input(Url(opt.input_name));

  // Use given option first, then use filter recommendation (it will probably
  // just recommend the value from the input platefile)

  if (opt.mode.empty())
    opt.mode = filter.mode(input.index_header().type());
  if (opt.tile_size == 0)
    opt.tile_size = filter.tile_size(input.default_tile_size());
  if (opt.filetype.empty())
    opt.filetype = filter.filetype(input.default_file_type());
  if (opt.pixel_format == VW_PIXEL_UNKNOWN)
    opt.pixel_format = filter.pixel_format(input.pixel_format());
  if (opt.channel_type == VW_CHANNEL_UNKNOWN)
    opt.channel_type = filter.channel_type(input.channel_type());

  PlateFile output(opt.output_name, opt.mode, opt.description, opt.tile_size, opt.filetype, opt.pixel_format, opt.channel_type);

  Transaction output_transaction_id = output.transaction_request("plate2plate, reporting for duty", -1);

  filter.init(output, input, input.transaction_cursor(), output_transaction_id);

  VW_ASSERT(input.num_levels() < 31, ArgumentErr() << "Can't handle plates deeper than 32 levels");

  int bottom_level = min(input.num_levels(), opt.bottom_level+1);

  for (int level = 0; level < bottom_level; ++level) {
    vw_out(InfoMessage) << "\nProcessing level " << level << " of " << bottom_level-1 << ".  ";
    TerminalProgressCallback tpc("plate.plate2plate.progress", "");
    vw::Timer timer( "\t    Processing time in seconds" );

    // The entire region contains 2^level tiles.
    int32 region_size = 1 << level;
    int subdivided_region_size = 1024;
    if (subdivided_region_size < region_size) subdivided_region_size = region_size;

    BBox2i full_region(0,0,region_size,region_size);
    std::list<BBox2i> boxes1 = bbox_tiles(full_region,
                                          subdivided_region_size,
                                          subdivided_region_size);

    vw_out(InfoMessage) << "Region"   << full_region << " has " << boxes1.size() << " bboxes\n";

    float region_counter = 0;
    BOOST_FOREACH( const BBox2i& region1, boxes1 ) {
      std::cout << "\n\t--> Sub-region: " << region1 << "\n";
      std::list<TileHeader> tiles;
      if ( !opt.skim_mode )
        tiles = input.search_by_region(level, region1, 0,
                                       input.transaction_cursor(), true);
      else
        tiles = input.search_by_region(level, region1, -1, -1, true);

      //      if (tiles.size() > 0)
      //      std::cout << "\t--> Region " << region1 << " has " << tiles.size() << " tiles.\n";
      std::ostringstream ostr;
      ostr << "\t    Converting " << tiles.size() << " tiles: ";
      tpc.set_progress_text(ostr.str());
      SubProgressCallback sub_progress(tpc,
                                       region_counter / boxes1.size(),
                                       (region_counter+1.0) / boxes1.size());
      BOOST_FOREACH( const TileHeader& tile, tiles ) {
        // ++n;
        // if (n % 100 == 0)
        //   std::cout << "n = " << n << "  --  "<< tile.col() << " " << tile.row() << " " << tile.level() << " " << tile.transaction_id() << "\n";
        filter(output, input, tile.col(), tile.row(), tile.level(),
               tile.transaction_id(), output_transaction_id);
        sub_progress.report_incremental_progress(1.0/tiles.size());
      }
      sub_progress.report_finished();
      region_counter++;
    }
    tpc.report_finished();
    output.sync();
  }

  filter.fini(output, input, output_transaction_id);

  output.transaction_complete(output_transaction_id, true);
}

// Blah blah boilerplate
int main(int argc, char *argv[])
{
  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    boost::to_lower(opt.filter);
    if (opt.filter == "identity") {
      Identity f;
      run(opt, f);
    } else if (opt.filter == "toast_dem") {
      ToastDem f;
      run(opt, f);
    }
  } catch (const Usage& e) {
    std::cout << e.what() << std::endl;
    return 1;
  } catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
