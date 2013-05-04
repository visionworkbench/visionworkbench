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


#include <vw/Core.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/FundamentalTypes.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/path.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main( int argc, char *argv[] ) {

  std::string blob_name;
  TransactionOrNeg transaction_id;

  po::options_description general_options("Writes out all inner images of a blob.\n\nGeneral Options");
  general_options.add_options()
    ("transaction_id,t", po::value(&transaction_id)->default_value(-1),
     "Transaction id to pull from blob.")
    ("help,h", "Display this help message.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("blob", po::value(&blob_name));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("blob", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <blob name>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    VW_OUT() << "An error occured while parsing command line arguments.\n\n";
    VW_OUT() << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    VW_OUT() << usage.str();
    return 1;
  }

  if( blob_name.empty() ) {
    std::cerr << "Error: must specify a platefile url.\n";
    VW_OUT() << usage.str();
    return 1;
  }

  // -------------------

  ReadBlob blob(blob_name);

  std::string blob_prefix = fs::path(blob_name).stem().string();

  VW_OUT() << "Searching for Transaction: " << transaction_id << "\n";
  TerminalProgressCallback tpc("plate.tools.blob_file_dump","Dumping:");
  tpc.report_progress(0);

  for ( ReadBlob::iterator rec_it = blob.begin(); rec_it != blob.end(); rec_it++ ) {
    // Get record and update the progress bar.
    BlobTileRecord rec = *rec_it;
    tpc.report_progress( double(rec_it.current_base_offset()) /
                         double(blob.size()) );

    // Check to see if this record has a TID that we want
    if (transaction_id >= 0 && rec.hdr.transaction_id() != transaction_id
        && !transaction_id.newest())
      continue;

    std::ostringstream ostr;
    ostr << blob_prefix << "_" << rec.hdr.transaction_id() << "_"
         << rec.hdr.level() << "_" << rec.hdr.row() << "_"
         << rec.hdr.col() << "." << rec.hdr.filetype();

    VW_OUT(VerboseDebugMessage,"plate.tools.blob_file_dump") << "Dumping: " << ostr.str() << "\n";

    // Here is where we actually write the data to file
    std::ofstream ofile(ostr.str().c_str(), std::ios::binary);
    VW_ASSERT(ofile.is_open(), IOErr() << "could not open dst file for writing (" << ostr.str() << ")");
    ofile.write(reinterpret_cast<char*>(&rec.data->operator[](0)), rec.data->size());
    VW_ASSERT(!ofile.fail(), IOErr() << ": failed to write to " << ostr.str());
    ofile.close();
  }
  tpc.report_finished();

  return 0;
}
