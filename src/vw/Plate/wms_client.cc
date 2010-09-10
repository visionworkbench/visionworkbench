// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//// Vision Workbench
#include <vw/Plate/WMSMessages.pb.h>
#include <vw/Plate/IndexService.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/common.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::platefile;
VW_DEFINE_EXCEPTION(Usage, Exception);

struct Options {
  float north, south, east, west; // Request area in lat lon degrees
  int32 width, height;
  int32 platefile_id;
};

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

struct RPC {
  boost::shared_ptr<AmqpRpcClient> client;
  boost::shared_ptr<WMSService>  service;
  RPC() {
    std::string queue_name = AmqpRpcClient::UniqueQueueName("wms_client");

    boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
    client.reset(  new AmqpRpcClient(conn, PLATE_EXCHANGE_NAMESPACE ".wms", queue_name, "wms") );
    service.reset( new WMSService::Stub(client.get() ) );
    client->bind_service(service, queue_name);
  }
};

VW_DEFINE_SINGLETON(rpc, RPC);

BBoxContainer adapt_bbox(const BBox2& in_box) {
  BBoxContainer bbox;
  bbox.set_origin_x(in_box.min().x());
  bbox.set_origin_y(in_box.min().y());
  bbox.set_width(in_box.width());
  bbox.set_height(in_box.height());
  return bbox;
}

void run(const Options& opt) {
  WMSTileRequest  q;
  WMSTileResponse a;

  q.set_platefile_id(opt.platefile_id);
  q.mutable_lonlat()->CopyFrom(adapt_bbox(BBox2(opt.west, opt.north, opt.east-opt.west, opt.south-opt.north)));
  q.mutable_pixels()->CopyFrom(adapt_bbox(BBox2(0, 0, opt.width, opt.height)));

  rpc_mutable().service->GetTile(rpc_mutable().client.get(), &q, &a, null_callback());

  vw_out(InfoMessage) << a.filename() << std::endl;
}

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Requests a wms tile");
  general_options.add_options()
    ("platefile,p", po::value(&opt.platefile_id),    "Request from this platefile.")
    ("north",       po::value(&opt.north),           "BBox (degrees)")
    ("south",       po::value(&opt.south),           "BBox (degrees)")
    ("west",        po::value(&opt.west),            "BBox (degrees)")
    ("east",        po::value(&opt.east),            "BBox (degrees)")
    ("width",       po::value(&opt.width), "Tile size (integer pixels)")
    ("height",      po::value(&opt.height), "Tile size (integer pixels)")
    ("help", "Display this help message");

  po::options_description options("Allowed Options");
  options.add(general_options);

  po::variables_map vm;

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-p <platefile_name>]" << std::endl;
  usage << general_options << std::endl;

  try {
    po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
    po::notify( vm );
  } catch (const po::error &e) {
    vw_throw(Usage() << "Error parsing input:\n\t" << e.what() << "\n" << usage.str());
  }

  if( vm.count("help") )
    vw_throw(Usage() << usage.str());

#define REQUIRE(x) do {\
  if( vm.count(x) != 1)\
    vw_throw(Usage() << usage.str() << "\nNeed a " x "\n");\
} while(0)

    REQUIRE("platefile");
    REQUIRE("north");
    REQUIRE("south");
    REQUIRE("east");
    REQUIRE("west");
    REQUIRE("width");
    REQUIRE("height");
}

int main(int argc, char** argv) {
  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    run(opt);
  } catch (const Usage& e) {
    std::cout << e.what() << std::endl;
    return 1;
  } catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}



