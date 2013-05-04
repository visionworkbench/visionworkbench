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


#include <iostream>
#include <fstream>

#include <vw/Core.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/detail/LocalIndex.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::platefile;

template <class T>
uint32 cast_uint32( T const& t ) {
  return boost::lexical_cast<uint32>(t);
}

struct Options {
  std::string plate_file_name;
  int32 blob_id;
};

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("Chops platefile into georeferenced squares.\n");
  general_options.add_options()
    ("help,h", "Display this help message");

  po::options_description positional("");
  positional.add_options()
    ("plate-file", po::value(&opt.plate_file_name))
    ("blob-id", po::value(&opt.blob_id)->default_value(-1));

  po::positional_options_description positional_desc;
  positional_desc.add("plate-file", 1);
  positional_desc.add("blob-id", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <local platefile dir> <blob id>\n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.plate_file_name.empty() || opt.blob_id < 0 )
    vw_throw( ArgumentErr() << "Missing input platefile or blob id!\n"
              << usage.str() << general_options );
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Work out how many files I have to process
    int32 count = 0;
    fs::recursive_directory_iterator end_itr;
    for ( fs::recursive_directory_iterator itr( fs::path(opt.plate_file_name + "/index" ) );
          itr != end_itr; ++itr ) {
      if ( fs::is_directory(itr->status()) )
        continue;
      count++;
    }

    // Structure to put results in. The key is the offset inside the
    // blob. The value is <level,row,col,tid>.
    typedef Vector<uint32,4> Vector4u;
    std::map<size_t,Vector4u> blob_contents;

    // Actually start opening the index and querying
    TerminalProgressCallback tpc("","");
    tpc.report_progress(0);
    double tpc_inc = 1.0 / double(count);
    for ( fs::recursive_directory_iterator itr( fs::path(opt.plate_file_name + "/index" ) );
          itr != end_itr; ++itr ) {
      if ( fs::is_directory(itr->status()) )
        continue;

      std::vector<std::string> split_vec;
      boost::split( split_vec, itr->path().string(),
                    boost::is_any_of("/") );
      size_t elements = split_vec.size();

      // Page size is normally 256 px by 256 px. I'm going to assume that.
      Vector<uint32,2> page_size(256,256);
      uint32 level = cast_uint32(split_vec[elements-3]);
      uint32 col   = cast_uint32(split_vec[elements-2]);
      uint32 row   = cast_uint32(split_vec[elements-1]);

      // Deserialize
      detail::LocalIndexPage page(itr->path().string(),
                                  level, col, row,
                                  page_size[0], page_size[1] );

      // Unfortunately we can't just iterator over the populated areas
      // of LocalIndexPage's sparsetable. Sparsetable doesn't keep the
      // keys per say. We could get the values (the populated Blob)
      // but we wouldn't know what tile location that corresponds to.

      for ( uint32 j = row; j < row+page_size[1]; j++ ) {
        for ( uint32 i = col; i < col+page_size[0]; i++ ) {
          try {
            BOOST_FOREACH( platefile::detail::IndexPage::value_type const& record,
                           page.multi_get( i, j, 0, platefile::detail::MAX_TRANSACTION ) ) {
              if ( record.second.blob_id() == opt.blob_id )
                blob_contents[record.second.blob_offset()] = Vector4u(level,j,i,record.first);
            }
          } catch ( const platefile::TileNotFoundErr& e) {}
        }
      }

      tpc.report_incremental_progress( tpc_inc );
    }
    tpc.report_finished();

    if ( blob_contents.empty() ) {
      VW_OUT() << "Index makes no mention of that blob ID\n";
      return 0;
    }

    // Print out the results of what was supposed to be in the blob.
    VW_OUT() << "Supposed contents of Blob = " << opt.blob_id << std::endl;
    VW_OUT() << "Offset | Level, Row, Col, Tid\n";
    for ( std::map<size_t,Vector4u>::const_iterator it = blob_contents.begin();
          it != blob_contents.end(); it++ ) {
      VW_OUT() << it->first << " " << it->second << std::endl;
    }

  } catch ( const ArgumentErr& e ) {
    VW_OUT() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "VW Error: " << e.what() << std::endl;
    return 1;
  } catch ( const std::bad_alloc& e ) {
    std::cerr << "Error: Ran out of Memory!" << std::endl;
    return 1;
  } catch ( const std::exception& e ) {
    std::cerr << "Error: " << e.what() <<  std::endl;
    return 1;
  }

  return 0;
}
