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


#include <vw/Plate/detail/Index.h>
#include <vw/Plate/detail/IndexPage.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/HTTPUtils.h>

#include <sstream>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::platefile;


VW_DEFINE_EXCEPTION(Usage, Exception);

struct Options {
  uint32 col, row, level;
  Url plate;
  bool verify;
  //std::string output;
};

boost::shared_ptr<ReadBlob> open_blob(const Url& url, uint32 blob_id) {
  static const boost::format blob_tmpl("%s/plate_%u.blob");
  static uint32 last_blob_id = std::numeric_limits<uint32>::max();
  static boost::shared_ptr<ReadBlob> last_blob;

  if (blob_id != last_blob_id) {
    boost::format blob_name(blob_tmpl);
    last_blob.reset(new ReadBlob(boost::str(blob_name % url.path() % blob_id)));
    last_blob_id = blob_id;
  }
  return last_blob;
}

void dump_tile(const Url& url, uint32 blob_id, uint64 blob_offset) {
  boost::shared_ptr<ReadBlob> blob = open_blob(url, blob_id);
  TileHeader hdr = blob->read_header(blob_offset);
  std::cerr << "   => TID=" << hdr.transaction_id() << " COL=" << hdr.col() << " ROW=" << hdr.row() << " LEVEL=" << hdr.level() << std::endl;
}

void run(const Options& opt) {
  boost::shared_ptr<detail::Index> index(detail::Index::construct_open(opt.plate));
  boost::shared_ptr<detail::IndexPage> page = index->page_request(opt.col, opt.row, opt.level);
  std::cout << "Loaded page at col=" << opt.col << " row=" << opt.row << " level=" << opt.level << std::endl
            << "Page contains " << page->sparse_size() << " entries." << std::endl;

    //typedef std::pair<uint32, IndexRecord> value_type;
    //typedef std::list<value_type> multi_value_type;
    //typedef google::sparsetable<multi_value_type>::nonempty_iterator nonempty_iterator;

  BOOST_FOREACH(const detail::IndexPage::multi_value_type& slot, std::make_pair(page->begin(), page->end())) {
    std::cout << "Loaded page slot with " << slot.size() << " entries" << std::endl;
    BOOST_FOREACH(const detail::IndexPage::value_type& elt, slot) {
      std::cout << "TID=" << elt.first << " BLOB=" << elt.second.blob_id() << " OFFSET=" << elt.second.blob_offset() << std::endl;
      if (opt.verify)
        dump_tile(opt.plate, elt.second.blob_id(), elt.second.blob_offset());
    }
  }
}

void process_args(Options& opt, int argc, char *argv[]) {
  po::options_description general_options("Dump out information about an index page.\n\nGeneral Options");
  general_options.add_options()
    ("col,c",         po::value(&opt.col),    "col to get page for")
    ("row,r",         po::value(&opt.row),    "row to get page for")
    ("level,l",       po::value(&opt.level),  "level to get page for")
    ("verify-blob",   po::bool_switch(&opt.verify)->default_value(false), "also look up blob entry")
    //("output-name,o", po::value(&opt.output), "Specify the base output directory")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("plate-file", po::value(&opt.plate));

  po::options_description options("");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("plate-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw(Usage() << "Could not parse command line arguments: " << usage.str());
  }

  if( vm.count("help") )
    vw_throw(Usage() << usage.str());

#define REQUIRE(arg)\
  if (!vm.count(#arg)) {\
    vw_throw(Usage() << "Error: Expected a " #arg " argument\n" << usage.str());\
  }\

  REQUIRE(plate-file);
  REQUIRE(col);
  REQUIRE(row);
  REQUIRE(level);

  if (opt.verify && opt.plate.scheme() != "file")
    vw_throw(Usage() << "Can only use --verify-blob with a file:// scheme");

  //if (opt.output.empty())
  //  opt.output = (opt.plate.scheme() == "file" ? opt.plate.path() : PlatefileUrl(opt.plate).name()) + ".tiles";
}

int main(int argc, char** argv) {
  Options opt;
  try {
    process_args(opt, argc, argv);
    run(opt);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
