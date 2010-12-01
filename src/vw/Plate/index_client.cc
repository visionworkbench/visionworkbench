// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//// Vision Workbench
#include <vw/Plate/IndexService.h>
#include <vw/Plate/Index.h>
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

  vw_out() << "Platefile: "
           << "ID["          << hdr.platefile_id()      << "] "
           << "Name["        << fs::path(url.path()).leaf()   << "] "
           << "Filename["    << index->platefile_name() << "] "
           << "Description[" << (hdr.has_description() ? hdr.description() : "No Description") << "] "
           << "MaxLevel["    << index->num_levels()-1   << "]"
           << std::endl;
}

void ListPlates(const Url& url) {

  IndexListRequest request;
  IndexListReply   reply;

  boost::scoped_ptr<IndexClient> client(conn(url));

  std::cerr << "Listing plates!" << std::endl;
  client->ListRequest(client.get(), &request, &reply, null_callback());

  vw_out() << "Got Plates:" << std::endl;
  std::copy(reply.platefile_names().begin(), reply.platefile_names().end(), std::ostream_iterator<std::string>(vw_out(), " "));
  vw_out() << std::endl;

  BOOST_FOREACH(const std::string& name, reply.platefile_names())
    PlateInfo(PlatefileUrl(url, name));
}

int main(int argc, char** argv) {
  Url url;

  po::options_description general_options("Runs a query against the index manager, or a specified platefile id");
  general_options.add_options()
    ("url,u", po::value(&url), "Run an info request against this platefile url.")
    ("help,h",  "Display this help message");

  po::options_description options("Allowed Options");
  options.add(general_options);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-p <platefile_name>]" << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if (fs::extension(url.path()) == ".plate") {
    // Run IndexInfoRequest
    PlateInfo(url);
  } else {
    // Run IndexListRequest
    ListPlates(url);
  }

  return 0;
}

