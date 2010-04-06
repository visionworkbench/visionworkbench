#include <vw/Plate/ToastDem.h>
#include <vw/Plate/PlateManager.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Core.h>
#include <boost/shared_ptr.hpp>
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

  inline void init(PlateFile& output, const PlateFile& input, int transaction_id) {
    return impl().init(output, input, transaction_id);
  }
  inline void fini(PlateFile& output, const PlateFile& input, int transaction_id) {
    return impl().fini(output, input, transaction_id);
  }

#   define lookup(name, type) type name(type data) const { return impl().name(data); }
    lookup(mode, string);
    lookup(tile_size, int);
    lookup(filetype, string);
    lookup(pixel_format, PixelFormatEnum);
    lookup(channel_type, ChannelTypeEnum);
#   undef lookup

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, int32 transaction_id) {
    impl()(output, input, col, row, level, transaction_id);
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

  inline void init(PlateFile& output, const PlateFile& /*input*/, int /*transaction_id*/) { output.write_request(); }
  inline void fini(PlateFile& output, const PlateFile& /*input*/, int /*transaction_id*/) { output.write_complete(); }

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, int32 transaction_id) {
    ImageView<PixelRGBA<double> > tile;
    TileHeader hdr = input.read(tile, col, row, level, transaction_id);
    output.write_update(tile, col, row, level, transaction_id);
  }
};

struct ToastDem : public FilterBase<ToastDem> {
  string mode(string) const { return "toast_dem"; }
  int tile_size(int)  const { return 32; }
  string filetype(string) const { return "toast_dem_v1"; }
  PixelFormatEnum pixel_format(PixelFormatEnum) const { return VW_PIXEL_SCALAR; }
  ChannelTypeEnum channel_type(ChannelTypeEnum) const { return VW_CHANNEL_UINT8; }

  inline void init(PlateFile& output, const PlateFile& input, int transaction_id) {
    output.write_request();

    // Write null tiles for the levels we don't have data for
    int level_difference = log(input.default_tile_size()/float(output.default_tile_size())) / log(2.) + 0.5;

    vw_out(InfoMessage, "plate.tools.plate2plate") << "Creating null tiles for a level difference of " << level_difference << std::endl;

    uint64 bytes;
    boost::shared_array<uint8> null_tile = toast_dem_null_tile(bytes);

    for (int level = 0; level < level_difference; ++level) {
      int region_size = 1 << level;
      for (int col = 0; col < region_size; ++col)
        for (int row = 0; row < region_size; ++row)
          output.write_update(null_tile, bytes, col, row, level, transaction_id);
    }
  }

  inline void fini(PlateFile& output, const PlateFile& /*input*/, int /*transaction_id*/) {
    output.write_complete();
  }

  struct DemWriter : public ToastDemWriter {
    PlateFile& platefile;
    DemWriter(PlateFile& output) : platefile(output) { }
    inline void operator()(const boost::shared_array<uint8> data, uint64 data_size, int32 dem_col, int32 dem_row, int32 dem_level, int32 transaction_id) const {
      platefile.write_update(data, data_size, dem_col, dem_row, dem_level, transaction_id);
    }
  };

  inline void operator()( PlateFile& output, const PlateFile& input, int32 col, int32 row, int32 level, int32 transaction_id) {
    DemWriter writer(output);
    make_toast_dem_tile(writer, input, col, row, level, transaction_id);
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

  string filter;

  Options() :
    tile_size(0), pixel_format(VW_PIXEL_UNKNOWN), channel_type(VW_CHANNEL_UNKNOWN) {}
};

VW_DEFINE_EXCEPTION(Usage, Exception);

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description options("Options");
  options.add_options()
    ("output-name,o",    po::value(&opt.output_name), "Specify the URL of the input platefile.")
    ("input-name,i",     po::value(&opt.input_name),  "Specify the URL of the output platefile.")
    ("file-type",        po::value(&opt.filetype),    "Output file type")
    ("mode",             po::value(&opt.mode),        "Output mode [toast, kml]")
    ("tile-size",        po::value(&opt.tile_size),   "Output size, in pixels")
    ("filter",           po::value(&opt.filter),      "Filters to run")
    ("help,h",           "Display this help message.");

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw(Usage() << "Error parsing input:\n\t"
             << e.what() << "\n" << options );
  }

  if (opt.output_name.empty() || opt.input_name.empty() || opt.filter.empty())
    vw_throw(Usage() << options);
}

template <typename FilterT>
void run(Options& opt, FilterBase<FilterT>& filter) {
  PlateFile input(opt.input_name);

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

  int transaction_id = output.transaction_request("plate2plate, reporting for duty", -1);

  filter.init(output, input, transaction_id);

  VW_ASSERT(input.num_levels() < 31, ArgumentErr() << "Can't handle plates deeper than 32 levels");

  for (int level = 0; level < input.num_levels(); ++level) {
    vw_out(InfoMessage) << "Processing level " << level << " of " << input.num_levels()-1 << std::endl;
    TerminalProgressCallback tpc("plate.plate2plate.progress", "");
    vw::Timer::Timer( "Processed in" );

    // The entire region contains 2^level tiles.
    int32 region_size = 1 << level;

    double step = 1./(region_size*region_size);
    tpc.print_progress();

    for (int32 i = 0; i < region_size; ++i) {
      for (int32 j = 0; j < region_size; ++j) {
        std::list<TileHeader> tiles;
        try {
          tiles = input.search_by_location(i, j, level, 0, std::numeric_limits<int>::max());
        } catch (const TileNotFoundErr&) { continue; }

        BOOST_FOREACH(const TileHeader& tile, tiles)
          filter(output, input, tile.col(), tile.row(), tile.level(), transaction_id);
        tpc.report_incremental_progress(step);
      }
    }
    tpc.report_finished();
    output.sync();
  }

  filter.fini(output, input, transaction_id);

  output.transaction_complete(transaction_id, true);
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
