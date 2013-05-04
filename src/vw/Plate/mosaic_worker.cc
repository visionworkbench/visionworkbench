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


#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>

// Protocols
#include <vw/Plate/Rpc.pb.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::platefile;




// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string queue_name;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value(""), "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("help,h", "Display this help message");

  po::options_description options("Allowed Options");
  options.add(general_options);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-q <queue name>]" <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  std::cout << "Opening AMQP connection.\n";
  AmqpConnection conn;

  conn.exchange_declare("ngt", "direct", true);
  conn.basic_publish("hello there", "ngt", "jobs");

  conn.queue_declare("test", true, true, true);
  conn.queue_bind("test", "ngt", "jobs");
  //  while(1)
    //    conn.basic_consume("test", "testtag", false);

  return 0;
}
