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


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/HTTPUtils.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  Url url;
  unsigned transaction_id;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("complete", po::value<unsigned>(&transaction_id), "Mark a transaction as complete, and update the read cursor..")
    ("sync", "Sync the platefile index to disk.")
    ("help,h", "Display this help message.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("url", po::value(&url));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("url", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <url>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("url") < 1 ) {
    std::cerr << "Error: must specify a platefile url.\n";
    std::cout << usage.str();
    return 1;
  }

  //------------------------- SET DEFAULT OPTIONS -----------------------------

  try {

    std::cout << "\nOpening plate file: " << url.string() << "\n";
    boost::shared_ptr<PlateFile> platefile =
      boost::shared_ptr<PlateFile>( new PlateFile(url) );

    if (vm.count("complete")) {
      platefile->transaction_resume(transaction_id);
      platefile->transaction_end(true);
    }

    if (vm.count("sync"))
      platefile->sync();

  }  catch (const vw::Exception& e) {
    std::cout << "An error occured: " << e.what() << "\nExiting.\n\n";
    return 1;
  }

  return 0;
}
