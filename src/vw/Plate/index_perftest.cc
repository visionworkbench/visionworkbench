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
    return 1;
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

  return 0;
}

