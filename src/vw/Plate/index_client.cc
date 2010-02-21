// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//// Vision Workbench
#include <vw/Plate/IndexService.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/common.h>

#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::platefile;

static void null_closure() {}

#if 0
const boost::shared_ptr<Blob> PlateModule::get_blob(const std::string& plate_filename, uint32 blob_id) const {
  std::ostringstream ostr;
  ostr << plate_filename << "/plate_" << blob_id << ".blob";
  const std::string& filename = ostr.str();

  BlobCache::const_iterator blob = blob_cache.find(filename);
  if (blob != blob_cache.end())
    return blob->second;

  boost::shared_ptr<Blob> ret( new Blob(filename, true) );
  blob_cache[filename] = ret;
  return ret;
}
#endif


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

std::string queue_name() {
  return AmqpRpcClient::UniqueQueueName("index_client");
}

struct RPC {
  boost::shared_ptr<AmqpRpcClient> client;
  boost::shared_ptr<IndexService>  service;
  RPC() {
    boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
    client.reset(  new AmqpRpcClient(conn, DEV_INDEX, queue_name(), "index") );
    service.reset( new IndexService::Stub(client.get() ) );
    client->bind_service(service, queue_name());
  }
};

VW_DEFINE_SINGLETON(rpc, RPC);

void PlateInfo(const std::string& name) {

  boost::shared_ptr<Index> index = Index::construct_open(std::string("pf://index/") + name);
  const IndexHeader& hdr = index->index_header();

  vw_out() << "Platefile: "
           << "ID["          << hdr.platefile_id()      << "] "
           << "Name["        << fs::path(name).leaf()   << "] "
           << "Filename["    << index->platefile_name() << "] "
           << "Description[" << (hdr.has_description() ? hdr.description() : "No Description") << "]"
           << std::endl;
}

void ListPlates() {

  IndexListRequest request;
  IndexListReply   reply;

  rpc_mutable().service->ListRequest(rpc_mutable().client.get(), &request, &reply, google::protobuf::NewCallback(&null_closure));

  vw_out() << "Got Plates:" << std::endl;
  std::copy(reply.platefile_names().begin(), reply.platefile_names().end(), std::ostream_iterator<std::string>(vw_out(), " "));
  vw_out() << std::endl;

  std::for_each(reply.platefile_names().begin(), reply.platefile_names().end(), boost::bind(&PlateInfo, _1));
}

int main(int argc, char** argv) {

  std::string name;

  po::options_description general_options("Runs a query against the index manager, or a specified platefile id");
  general_options.add_options()
    ("platefile,p", po::value(&name), "Run an info request against this platefile id.")
    ("help", "Display this help message");

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

  if (vm.count("platefile")) {
    // Run IndexInfoRequest
    PlateInfo(name);
  } else {
    // Run IndexListRequest
    ListPlates();
  }

  return 0;
}

