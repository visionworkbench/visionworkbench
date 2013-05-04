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
#include <vw/Core/Stopwatch.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Log.h>
#include <signal.h>
#include <boost/format.hpp>

#include <google/protobuf/descriptor.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

// ------------------------------ SIGNAL HANDLER -------------------------------

volatile bool process_messages = true;
volatile bool force_sync = false;

void sig_unexpected_shutdown(int sig_num) {
  signal(sig_num, SIG_IGN);
  process_messages = false;
  signal(sig_num, sig_unexpected_shutdown);
}

void sig_sync(int sig_num) {
  signal(sig_num, SIG_IGN);
  force_sync = true;
  signal(sig_num, sig_sync);
}

struct Options {
  Url url;
  std::string root;
  float sync_interval;
  bool debug;
  bool help;
};

VW_DEFINE_EXCEPTION(Usage, Exception);

void process_args(Options& opt, int argc, char *argv[]) {
  po::options_description general_options("Runs a master index manager.\n\nGeneral Options:");
  general_options.add_options()
    ("url",             po::value(&opt.url),                               "Url to listen on")
    ("debug",           po::bool_switch(&opt.debug)->default_value(false), "Allow server to die.")
    ("help,h",          po::bool_switch(&opt.help)->default_value(false),  "Display this help message")
    ("sync-interval,s", po::value(&opt.sync_interval)->default_value(60.),
     "Specify the time interval (in minutes) for automatically synchronizing the index to disk.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("root-directory", po::value(&opt.root));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("root-directory", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " --url <url> root_directory" << std::endl << std::endl;
  usage << general_options << std::endl;

  if(opt.help)
    vw_throw(Usage() << usage.str());

  if( vm.count("root-directory") != 1 ) {
    vw_throw(Usage() << usage.str()
                     << "\n\nError: must specify a root directory that contains plate files!");
  }

  if ( vm.count("url") != 1 ) {
    vw_throw(Usage() << usage.str()
                     << "\n\nMust specify a url to listen on");
  }
}

int main(int argc, char** argv) {
  Options opt;
  try {
    process_args(opt, argc, argv);
  } catch (const Usage& u) {
    std::cerr << u.what() << std::endl;
    ::exit(EXIT_FAILURE);
  }

  // Install Unix Signal Handlers.  These will help us to gracefully
  // recover and salvage the index under most unexpected error
  // conditions.
  signal(SIGINT,  sig_unexpected_shutdown);
  signal(SIGTERM,  sig_unexpected_shutdown);
  signal(SIGUSR1, sig_sync);

  // Start the server task in another thread
  RpcServer<IndexServiceImpl> server(opt.url, new IndexServiceImpl(opt.root));
  server.set_debug(opt.debug);

  VW_OUT(InfoMessage) << "Starting index server\n\n";
  uint64 sync_interval_us = uint64(opt.sync_interval * 60000000);
  uint64 t0 = Stopwatch::microtime(), t1;
  uint64 next_sync = t0 + sync_interval_us;

  size_t win = 0, lose = 0, draw = 0, total = 0;
  boost::format status("qps[%7.1f]   total[%9u]   server_err[%9u]   client_err[%9u]\r");

  while(process_messages) {
    if (server.error()) {
      VW_OUT(InfoMessage) << "IO Thread has terminated with message: " << server.error() << "\n";
      break;
    }

    bool should_sync = force_sync || (Stopwatch::microtime() >= next_sync);

    if (should_sync) {
      VW_OUT(InfoMessage) << "\nStarting sync to disk. (" << (force_sync ? "manual" : "auto") << ")\n";
      uint64 s0 = Stopwatch::microtime();
      server.impl()->sync();
      uint64 s1 = Stopwatch::microtime();
      next_sync = s1 + sync_interval_us;
      VW_OUT(InfoMessage) << "Sync complete (took " << float(s1-s0) / 1e6  << " seconds).\n";
      force_sync = false;
    }

    t1 = Stopwatch::microtime();

    size_t win_dt, lose_dt, draw_dt, total_dt;
    {
      ThreadMap::Locked stats = server.stats();
      win_dt = stats.get("msgs");
      lose_dt  = stats.get("server_error");
      draw_dt  = stats.get("client_error");
      stats.clear();
    }
    total_dt = win_dt + lose_dt + draw_dt;
    win   += win_dt;
    lose  += lose_dt;
    draw  += draw_dt;
    total += total_dt;

    float dt = float(t1 - t0) / 1e6f;
    t0 = t1;

    VW_OUT(InfoMessage)
      << status % (float(total_dt)/dt) % total % lose % draw
      << std::flush;

    Thread::sleep_ms(500);
  }

  VW_OUT(InfoMessage) << "\nShutting down the index service safely.\n";
  server.stop();
  server.impl()->sync();

  return 0;
}

