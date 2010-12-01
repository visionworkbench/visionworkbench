// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// tiles2plate.cc
///
/// Converts a nested set of directories and tiles to a plate file.

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Plate/PlateFile.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct Tile {
  fs::path m_path;
  int32 m_col, m_row, m_level;

  Tile(const fs::path &path, int32 col = -1, int32 row = -1, int32 level = -1)
    : m_path(path), m_col(col), m_row(row), m_level(level)
  { }

  inline bool valid() {
    return ( (m_row != -1 && m_col != -1 && m_level != -1) &&
             fs::exists(m_path) );
  }
};

class TilePathDecoder {
public:
  TilePathDecoder() {}
  virtual ~TilePathDecoder() {}
  virtual Tile decode( const fs::path &path ) = 0;
};

class ToastDecoder : public TilePathDecoder {
  fs::path m_base_path;

public:
  ToastDecoder( fs::path & base_path )
    : m_base_path(base_path) { }
  virtual ~ToastDecoder() { }

  virtual Tile decode( const fs::path &full_path ) {
    fs::path partial_path(full_path.string().substr(m_base_path.string().length() + 1));
    Tile t(full_path);

    partial_path.replace_extension();
    //std::cout << partial_path << std::endl;
    for( fs::path::iterator itr = partial_path.begin();
         itr != partial_path.end();
         ++itr ) {
      std::istringstream in(*itr);
      t.m_level = t.m_col;
      t.m_col = t.m_row;
      if ( !(in >> t.m_row) )
        t.m_row = -1;
    }

    return t;
  }
};

class GigapanDecoder : public TilePathDecoder {
public:
  GigapanDecoder() {}
  virtual ~GigapanDecoder() {}

  virtual Tile decode( const fs::path &full_path ) {
    Tile t(full_path);
    if( full_path.extension() == ".info" )
      return t;

    std::string tile_name = fs::path(full_path).replace_extension().leaf();

    if( tile_name[0] != 'r' )
      return t;

    tile_name.erase(0, 1);

    t.m_level = boost::numeric_cast<int32>(tile_name.length());
    if(t.m_level == 0) {
      t.m_row = t.m_col = 0;
      return t;
    }

    int32 col = 0, row = 0;
    for( std::string::iterator itr = tile_name.begin();
         itr != tile_name.end();
         ++itr ) {
      row <<= 1; col <<= 1;

      switch( *itr ) {
      case '0':
        break;
      case '1':
        col |= 1;
        break;
      case '2':
        row |= 1;
        break;
      case '3':
        row |= 1;
        col |= 1;
        break;
      default:
        return t;
      }
    }

    t.m_col = col;
    t.m_row = row;
    return t;
  }
};

class TileHandler {
public:
  TileHandler() {}
  virtual ~TileHandler() {}
  virtual bool operator() ( const Tile & tile ) = 0;
};

class FirstHandler : public TileHandler {
  fs::path m_path;
  bool m_have_path;

public:
  FirstHandler() : m_have_path(false) {}
  virtual ~FirstHandler() {}
  virtual bool operator() ( const Tile & tile ) {
    m_path = tile.m_path;
    m_have_path = true;
    return false;
  }
  bool have_path() { return m_have_path; }
  fs::path get_path() { return m_path; }
};

template <class ViewT>
class PlateHandler : public TileHandler {
  boost::shared_ptr<PlateFile> m_platefile;
  Transaction m_write_transaction_id;

public:
  PlateHandler( boost::shared_ptr<PlateFile> &platefile, Transaction write_transaction_id )
    : m_platefile(platefile), m_write_transaction_id(write_transaction_id)
  { }
  virtual ~PlateHandler() {}

  virtual bool operator() ( const Tile &tile ) {
    ViewT image(tile.m_path.string());
    std::cout << "\t--> Writing tile [ "
              << tile.m_col << " " << tile.m_row << " " << tile.m_level << " ] "
              << " (" << tile.m_path.leaf() << ")"
              << std::endl;
    m_platefile->write_request();
    m_platefile->write_update(image,
                              tile.m_col, tile.m_row,
                              tile.m_level, m_write_transaction_id);
    m_platefile->write_complete();
    return true;
  }
};

template <class TileHandlerT>
void iterate_over_tiles( fs::path &dir_path,
                         boost::shared_ptr<TilePathDecoder> &decoder,
                         TileHandlerT &handler) {

  fs::recursive_directory_iterator end_itr;
  for( fs::recursive_directory_iterator itr( dir_path );
       itr != end_itr;
       ++itr )
  {
    if( fs::is_directory(itr->status()) )
      continue;

    Tile t = decoder->decode( itr->path() );
    if( t.valid() ) {
      if( !handler(t) ) return;
    }
  }
}

int main( int argc, char *argv[] ) {
  std::string tile_format;
  std::string output_file_name;
  std::string output_file_type;
  std::string tile_directory_name;

  po::options_description general_options("Turns tile trees into a plate file.\n\nGeneral Options");
  general_options.add_options()
    ("format,f", po::value<std::string>(&tile_format)->default_value("toast"), "Tile directory format: toast, gigapan")
    ("output-name,o", po::value<std::string>(&output_file_name), "Specify the output plate file name.")
    ("file-type", po::value<std::string>(&output_file_type), "Output file type (png is used by default)")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("tile-directory", po::value<std::string>(&tile_directory_name));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("tile-directory", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <tile_directory>" << std::endl << std::endl
        << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch( po::error &e) {
    std::cout << "An error occurred while parsing command line options.\n\n";
    std::cout << usage.str();
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( tile_format != "toast" && tile_format != "gigapan" ) {
    vw_throw(ArgumentErr() << "Unknown fromat passed in using --format: " << tile_format);
  }

  if( vm.count("tile-directory") != 1 ) {
    std::cerr << "Error: must specify an input tile directory!" << std::endl << std::endl
              << usage.str();
    return 1;
  }

  if( !fs::exists(tile_directory_name) ) {
    std::cerr << "Error: input tile directory does not exist!" << std::endl << std::endl
              << usage.str();
    return 1;
  }

  fs::path tile_directory( tile_directory_name );

  if( !vm.count("output-name") ) {
    output_file_name = fs::path( tile_directory ).replace_extension("plate").string();
  }

  std::cout << "Reading tiles in '" << tile_format << "' format"
            << " from " << tile_directory_name
            << " into " << output_file_name << std::endl;

  boost::shared_ptr<TilePathDecoder> decoder;
  if ( tile_format == "toast" )
    decoder = boost::shared_ptr<TilePathDecoder>( new ToastDecoder(tile_directory) );
  else {
    decoder = boost::shared_ptr<TilePathDecoder>( new GigapanDecoder() );
  }

  FirstHandler first_handler;
  iterate_over_tiles( tile_directory, decoder, first_handler );

  if ( !first_handler.have_path() ) {
    std::cerr << "Error: no tiles in tile directory matching " << tile_format << " tile scheme!"
             << std::endl << std::endl
             << usage.str();
    return 1;
  }

  DiskImageResource *rsrc = DiskImageResource::open(first_handler.get_path().string());
  PixelFormatEnum pixel_format = rsrc->pixel_format();
  ChannelTypeEnum channel_type = rsrc->channel_type();
  int32 tile_size = rsrc->rows();

  if (pixel_format == VW_PIXEL_GRAY)
    pixel_format = VW_PIXEL_GRAYA;
  if (pixel_format == VW_PIXEL_RGB)
    pixel_format = VW_PIXEL_RGBA;
  delete rsrc;

  if( !vm.count("file-type") ) {
    if (channel_type == VW_CHANNEL_FLOAT32)
      output_file_type = "tif";
    else
      output_file_type = "png";
  }

  try {
    std::cout << "Opening plate file: " << output_file_name << std::endl;
    boost::shared_ptr<PlateFile> platefile =
      boost::shared_ptr<PlateFile>( new PlateFile(output_file_name, "gigapan", "",
                                                  tile_size, output_file_type,
                                                  pixel_format, channel_type) );

    std::vector<TileHeader> empty_tile_list;
    Transaction write_transaction_id =
      platefile->transaction_request("Writing tiles from tile tree " + tile_directory_name, -1);

    switch(pixel_format) {
    case VW_PIXEL_GRAY:
      switch(channel_type) {
      case VW_CHANNEL_UINT8: {
          PlateHandler<DiskImageView<PixelGray<uint8> > > plate_handler( platefile, write_transaction_id );
          iterate_over_tiles( tile_directory, decoder, plate_handler );
        }
        break;
      case VW_CHANNEL_FLOAT32: {
          PlateHandler<DiskImageView<PixelGray<float> > > plate_handler( platefile, write_transaction_id );
          iterate_over_tiles( tile_directory, decoder, plate_handler );
        }
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by tiles2plate.\n");
      }
      break;

    case VW_PIXEL_GRAYA:
      switch(channel_type) {
      case VW_CHANNEL_UINT8: {
          PlateHandler<DiskImageView<PixelGrayA<uint8> > > plate_handler( platefile, write_transaction_id );
          iterate_over_tiles( tile_directory, decoder, plate_handler );
        }
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by tiles2plate.\n");
      }
      break;

    case VW_PIXEL_RGB:
      switch(channel_type) {
      case VW_CHANNEL_UINT8: {
          PlateHandler<DiskImageView<PixelRGB<uint8> > > plate_handler( platefile, write_transaction_id );
          iterate_over_tiles( tile_directory, decoder, plate_handler );
        }
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by tiles2plate.\n");
      }
      break;

    case VW_PIXEL_RGBA:
      switch(channel_type) {
      case VW_CHANNEL_UINT8: {
          PlateHandler<DiskImageView<PixelRGBA<uint8> > > plate_handler( platefile, write_transaction_id );
          iterate_over_tiles( tile_directory, decoder, plate_handler );
        }
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by tiles2plate.\n");
      }
      break;

    default:
      vw_throw(ArgumentErr() << "Platefile contains a pixel type not supported by tiles2plate.\n");
    }

    // update_read_cursor == true below.
    platefile->transaction_complete( write_transaction_id, true );
  } catch (vw::Exception &e) {
    std::cout << "An error occured: " << e.what() << "\nExiting\n\n";
  }

  return 0;
}
