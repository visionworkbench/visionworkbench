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


//// Vision Workbench
#include <vw/Plate/IndexService.h>
#include <vw/Plate/detail/Index.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/Rpc.h>
#include <vw/Core/Log.h>

#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::platefile;
using namespace vw::platefile::detail;

#define VW_DEFINE_SINGLETON(name, klass) \
  namespace { \
    vw::RunOnce name ## _once = VW_RUNONCE_INIT; \
    boost::shared_ptr<klass> name ## _ptr; \
    void init_ ## name() { \
      name ## _ptr = boost::shared_ptr<klass>(new klass()); \
    } \
    void kill_ ## name() { \
      init_ ## name(); \
    } \
  } \
  const klass& name() { \
    name ## _once.run( init_ ## name ); \
    return *name ## _ptr; \
  } \
  klass& name ## _mutable() { \
    name ## _once.run( init_ ## name ); \
    return *name ## _ptr; \
  }

typedef RpcClient<IndexService> IndexClient;

IndexClient* conn(const Url& url) {
  std::cerr << "Connecting to index server at " << url.string() << std::endl;
  return new IndexClient(url);
}

void PlateInfo(const Url& url) {

  boost::shared_ptr<Index> index = Index::construct_open(url);
  const IndexHeader& hdr = index->index_header();

  VW_OUT() << "Platefile:\n"
           << "\tID["          << hdr.platefile_id()      << "]\n"
           << "\tName["        << fs::path(url.path()).filename()   << "]\n"
           << "\tFilename["    << index->platefile_name() << "]\n"
           << "\tDescription[" << (hdr.has_description() ? hdr.description() : "No Description") << "]\n"
           << "\tMaxLevel["    << index->num_levels()-1   << "]\n"
           << "\tTransactionRead[" << hdr.transaction_read_cursor() << "]\n"
           << "\tTransactionWrite[" << hdr.transaction_write_cursor() << "]\n"
           << "\tType[" << hdr.type() << "]"
           << std::endl;
}

void ListPlates(const Url& url) {

  IndexListRequest request;
  IndexListReply   reply;

  boost::scoped_ptr<IndexClient> client(conn(url));

  std::cerr << "Listing plates!" << std::endl;
  client->ListRequest(client.get(), &request, &reply, null_callback());

  VW_OUT() << "Got Plates:" << std::endl;
  std::copy(reply.platefile_names().begin(), reply.platefile_names().end(), std::ostream_iterator<std::string>(::vw::vw_out(), " "));
  VW_OUT() << std::endl;

  BOOST_FOREACH(const std::string& name, reply.platefile_names())
    PlateInfo(PlatefileUrl(url, name));
}

int main(int argc, char** argv) {
  Url url;

  po::options_description general_options("Runs a query against the index manager, or a specified platefile id");
  general_options.add_options()
    ("url,u", po::value(&url), "Run an info request against this platefile url.")
    ("help,h",  "Display this help message");

  po::positional_options_description p;
  p.add("url", 1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(general_options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <url> " << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  } else if ( url == Url() ) {
    std::cout << "URL is missing.\n\n" << usage.str();
    return 0;
  }

  if (fs::path(url.path()).extension().string() == ".plate") {
    // Run IndexInfoRequest
    PlateInfo(url);
  } else {
    // Run IndexListRequest
    ListPlates(url);
  }

  return 0;
}

