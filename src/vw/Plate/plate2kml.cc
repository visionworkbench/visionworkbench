// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <kml/dom.h>
#include <kml/engine/kml_file.h>
#include <kml/base/file.h>

struct Options {
  std::string output_name, url_name;
};

kmldom::LatLonAltBoxPtr
create_latlonaltbox( double south, double west, double degree_size,
                     kmldom::KmlFactory* factory ) {
  kmldom::LatLonAltBoxPtr latlonaltbox =
    factory->CreateLatLonAltBox();
  latlonaltbox->set_north( south + degree_size );
  latlonaltbox->set_south( south );
  latlonaltbox->set_east( west + degree_size );
  latlonaltbox->set_west( west );
  return latlonaltbox;
}

kmldom::LatLonBoxPtr
create_latlonbox( double south, double west, double degree_size,
                  kmldom::KmlFactory* factory ) {
  kmldom::LatLonBoxPtr latlonbox =
    factory->CreateLatLonBox();
  latlonbox->set_north( south + degree_size );
  latlonbox->set_south( south );
  latlonbox->set_east( west + degree_size );
  latlonbox->set_west( west );
  return latlonbox;
}

kmldom::LodPtr
create_lod( ssize_t min, ssize_t max, kmldom::KmlFactory* factory ) {
  kmldom::LodPtr lod = factory->CreateLod();
  lod->set_minlodpixels( min );
  lod->set_maxlodpixels( max );
  return lod;
}

std::string prefix_from_location( ssize_t col, ssize_t row, ssize_t level ) {
  std::string result = "";
  size_t count = 0;
  for ( ssize_t i = level; i >= 0; i-- ) {
    if ( (count % 3 == 0) && ( count != 0 ) )
      result += "/";
    size_t mask = 0x01 << i;
    result += ( ( ( col & mask  ) >> i ) |
                ( ( ( row & mask ) >> i ) << 1 ) ) + 48;
    count++;
  }
  return result;
}

void draw_kml_level ( std::string const& base_folder, ssize_t level,
                      BBox2i const& region, TransactionOrNeg id,
                      boost::shared_ptr<PlateFile> plate ) {
  // Do not delete this pointer. LibKML is managing this itself somewhere.
  kmldom::KmlFactory* factory = kmldom::KmlFactory::GetFactory();

  // Searching for transactions
  std::list<TileHeader> tile_records;
  tile_records = plate->search_by_region(level,region,id,id,0,false);

  // Drawing Level
  BOOST_FOREACH( TileHeader const& t, tile_records ) {
    fs::path path( base_folder + "/" +
                   prefix_from_location( t.col(), t.row(), t.level() ) );
    fs::create_directory( path.parent_path() );
    std::string path_prefix = path.string();
    plate->read_to_file(path_prefix, t.col(), t.row(),
                        t.level(), id, true);

    // Deciding math for lat lon range
    double deg_delta = 360.0 / double( 1 << t.level() );
    double lon_west = -180.0 + t.col()*deg_delta;
    double lat_south = 180 - deg_delta*(t.row()+1);

    std::list<TileHeader> links =
      plate->search_by_region(level+1,region*2,id,id,0,false);

    { // Producing this layer's kml
      kmldom::FolderPtr folder = factory->CreateFolder();
      folder->set_name( path.string() );

      // Create Network Links
      BOOST_FOREACH( TileHeader const& lower, links ) {
        double ldeg_delta = 360.0 / double( 1 << lower.level() );
        double llon_west = -180.0 + lower.col()*ldeg_delta;
        double llat_south = 180 - ldeg_delta*(lower.row()+1);

        kmldom::NetworkLinkPtr netlink = factory->CreateNetworkLink();
        kmldom::RegionPtr region = factory->CreateRegion();
        region->set_lod( create_lod( 128, -1, factory ) );
        region->set_latlonaltbox( create_latlonaltbox( llat_south, llon_west,
                                                       ldeg_delta, factory ) );
        netlink->set_region(region);
        fs::path link_path( base_folder + "/" +
                            prefix_from_location( lower.col(), lower.row(),
                                                  lower.level() ) );
        kmldom::LinkPtr linkptr = factory->CreateLink();
        linkptr->set_viewrefreshmode(kmldom::VIEWREFRESHMODE_ONREGION);
        if ( ( lower.level() % 3 == 0 ) && lower.level() != 0 ) {
          // Pulling out the directory change
          std::string link_string = path.filename()+"/"+link_path.filename()+".kml";
          linkptr->set_href( link_string );
        } else {
          linkptr->set_href( link_path.replace_extension(".kml").filename() );
        }
        netlink->set_link(linkptr);
        folder->add_feature( netlink );
      }

      // Create ground overlay
      kmldom::GroundOverlayPtr goverlay = factory->CreateGroundOverlay();
      kmldom::RegionPtr region = factory->CreateRegion();
      if ( level == 1 ) {
        // Top layerish. Layer 0 is never drawn because it's
        // boundaries are illegal (Lat range too large). So this step
        // is to make sure layer 1 can always be seen when zoomed out.
        region->set_lod( create_lod( 1, 513, factory ) );
      } else if ( links.size() == 0 ) {
        // End of branch
        region->set_lod( create_lod( 128, -1, factory ) );
      } else {
        // Original code went to 1024
        region->set_lod( create_lod( 128, 513, factory ) );
      }
      region->set_latlonaltbox( create_latlonaltbox( lat_south, lon_west,
                                                     deg_delta, factory ) );
      goverlay->set_region( region );
      goverlay->set_latlonbox( create_latlonbox( lat_south, lon_west,
                                                 deg_delta, factory ) );
      kmldom::IconPtr icon = factory->CreateIcon();
      icon->set_href(path.replace_extension(plate->default_file_type()).filename());
      goverlay->set_icon( icon );
      goverlay->set_draworder( 50 );
      folder->add_feature( goverlay );

      // Writing
      kmldom::KmlPtr kml = factory->CreateKml();
      kml->set_feature( folder );
      kmlengine::KmlFilePtr kml_file = kmlengine::KmlFile::CreateFromImport(kml);
      std::string xml_data;
      kml_file->SerializeToString(&xml_data);
      kmlbase::File::WriteStringToFile(xml_data,path.replace_extension(".kml").string());
    } // end kml writing

    // Recursive call
    BOOST_FOREACH( TileHeader const& lower, links ) {
      draw_kml_level( base_folder, lower.level(),
                      BBox2i( lower.col(), lower.row(), 1, 1 ),
                      id, plate );
    }
  }
}

void handle_arguments(int argc, char* argv[], Options& opt) {
  po::options_description general_options("Extract plate into KML tiles.");
  general_options.add_options()
    ("output_name,o", po::value(&opt.output_name),
     "Output name for the KML result.")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("url", po::value(&opt.url_name), "");

  po::options_description options("");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("url",-1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate_filename> [options]\n";

  if ( vm.count("help") || opt.url_name.empty() )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( *(opt.url_name.rbegin()) == '/' )
    opt.url_name = opt.url_name.substr(0,opt.url_name.size()-1);
  if ( opt.output_name.empty() ) {
    size_t folder_divide, ext_divide;
    folder_divide = opt.url_name.rfind("/");
    ext_divide = opt.url_name.rfind(".");
    opt.output_name =
      opt.url_name.substr(folder_divide+1,ext_divide-folder_divide-1);
  }
}

int main( int argc, char *argv[] ) {
  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    std::cout << "Opening: " << opt.url_name << "\n";
    boost::shared_ptr<PlateFile> platefile( new PlateFile(opt.url_name) );

    // Create out output directory and chdir
    fs::create_directory( opt.output_name );

    // Recursive call that goes depth first in generation.
    draw_kml_level( opt.output_name, 0, BBox2i(0,0,1,1),
                    TransactionOrNeg(-1), platefile);

    // Rename 0.kml to be the opening kml
    fs::rename( fs::path(opt.output_name) / "0.kml",
                fs::path(opt.output_name) / (opt.output_name+".kml") );

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
