// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
using namespace vw;
using namespace vw::platefile;

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
    } catch (po::error &e) {
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

    BBox2i region(opt.minx, opt.miny, opt.width, opt.height);
    TransactionOrNeg id(-1);
    std::list<TileHeader> tile_records;
    tile_records = plate->search_by_region(opt.level,region,id,id,0,false);

    printf("[\n");
    std::list<TileHeader>::iterator i;
    for (i=tile_records.begin(); i != tile_records.end(); i++) {
        if (i != tile_records.begin()) { std::cout << ",\n"; }
        TileHeader t = *i;
        printf("{'level':%u,'col':%u,'row':%u}", t.level(), t.col(), t.row() );
    }
    printf("\n]\n");

    return 0;
}
