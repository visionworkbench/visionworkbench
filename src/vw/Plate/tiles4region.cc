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
#include <vw/Plate/TileManipulation.h>
#include <vw/Plate/PlateManager.h>
#include <vw/Cartography/GeoReference.h>
using namespace vw;
using namespace vw::platefile;
using namespace vw::cartography;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
// #include <boost/foreach.hpp>

struct Options {
  std::string plate_url;
  uint64 minx, miny, width, height, level;
};

void handle_arguments(int argc, char* argv[], Options& opt) {
  po::options_description general_options("Get the tiles from a platefile at a particular level and region");
  general_options.add_options()
    ("level,l", po::value(&opt.level), "Level of the plate to search at")
    //("transaction,t")
    ("help,h", "Get help");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("plate_url", po::value(&opt.plate_url), "")
    ("minx", po::value(&opt.minx), "")
    ("miny", po::value(&opt.miny), "")
    ("width", po::value(&opt.width), "")
    ("height", po::value(&opt.height), "");

  po::options_description options("");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("plate_url", 1);
  p.add("minx", 1);
  p.add("miny", 1);
  p.add("width", 1);
  p.add("height", 1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate url> minx miny width height\n";

  if ( vm.count("help")  ||
       ! vm.count("level") || ! vm.count("plate_url") ||
       ! vm.count("minx") || ! vm.count("miny") || ! vm.count("width") || ! vm.count("height")
       ){
    vw_throw( ArgumentErr() << usage.str() << general_options );
  }
}


int main( int argc, char *argv[]) {
  Options opt;
  try {
    handle_arguments( argc, argv, opt);
  } catch (const Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  boost::shared_ptr<PlateFile> plate( new PlateFile(opt.plate_url) );


  if (opt.level >= plate->index_header().num_levels()) {
    std::cout << "{\n";
    std::cout << "  \"ok\": false,";
    std::cout << "  \"error\": \"TOO DEEP\"\n";
    std::cout << "}";
    return 1;
  } else {

    PlateManager<vw::PixelGrayA<uint8> >* pm =
      PlateManager<vw::PixelGrayA<uint8> >::make( plate->index_header().type(), plate);
    vw::cartography::GeoReference georef = pm->georeference(opt.level);

    BBox2i region(opt.minx, opt.miny, opt.width, opt.height);
    TransactionOrNeg id(-1);
    std::list<TileHeader> tile_records;
    tile_records = plate->search_by_region(opt.level,region,TransactionRange(id));

    printf("{ \"ok\": true,\n");
    printf("\"result\": [\n");
    std::list<TileHeader>::iterator i;
    for (i=tile_records.begin(); i != tile_records.end(); i++) {
      if (i != tile_records.begin()) { std::cout << ",\n"; }

      const TileHeader& t = *i;

      int tile_size = plate->default_tile_size();
      Vector2 corner1 = georef.point_to_lonlat(georef.pixel_to_point( Vector2(t.col() * tile_size, t.row() * tile_size ) - Vector2(0.5,0.5) ));
      Vector2 corner2 = georef.point_to_lonlat(georef.pixel_to_point( Vector2((t.col()+1) * tile_size - 1, (t.row()+1) * tile_size - 1 ) + Vector2(0.5,0.5) ));
      float west = corner1[0];
      float east = corner2[0];
      float north = corner1[1];
      float south = corner2[1];

      printf("{\"level\":%u,\"col\":%u,\"row\":%u,\"west\":%f, \"east\":%f, \"north\":%f, \"south\":%f, \"filetype\":\"%s\"}",
             t.level(), t.col(), t.row(), west, east, north, south, t.filetype().c_str() );
    }
    printf("\n]\n");
    printf("}");
  }

  return 0;
}
