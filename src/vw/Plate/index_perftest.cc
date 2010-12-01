// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Rpc.h>
#include <vw/Plate/IndexService.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Stopwatch.h>
#include  <iostream>

#include <boost/scoped_ptr.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::platefile;

const static size_t BATCH_SIZE = 5000;

typedef RpcClient<IndexService> IndexClient;

IndexClient* conn(const Url& url) {
  std::cerr << "Connecting to index server at " << url.string() << std::endl;
  return new IndexClient(url);
}

int main(int argc, char** argv) {
  Url url;

  po::options_description general_options("AMQP Performance Test Program");
  general_options.add_options()
    ("url,u", po::value(&url), "Run requests against this index url.")
    ("help,h", "Display this help message");

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(general_options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << "\n\n";
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  boost::scoped_ptr<IndexClient> client(conn(url));

  uint64 t0, t1;

  while (1) {
    t0 = Stopwatch::microtime();

    for (size_t i = 0; i < BATCH_SIZE; i++) {
      IndexTestRequest request;
      request.set_value(i);

      IndexTestReply response;
      client->TestRequest(client.get(), &request, &response, null_callback());

      if (i != response.value())
        std::cerr << "Error: IndexTestMessage failed!\n";
    }

    t1 = Stopwatch::microtime();
    std::cout << 1000000. * float(BATCH_SIZE) / float(t1-t0) << " msg/s" << std::endl;
  }
}

