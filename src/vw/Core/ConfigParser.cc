#include <iostream>
#include <fstream>
#include <istream>

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;


#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Settings.h>
#include <vw/Core/ConfigParser.h>

using boost::shared_ptr;
using namespace std;
using namespace vw;

namespace {
  int32 name2level(string name) {
    typedef map<string, MessageLevel> MapType;
    static MapType name_map;

    if (name_map.empty()) {
      name_map["InfoMessage"]         = InfoMessage;
      name_map["ErrorMessage"]        = ErrorMessage;
      name_map["WarningMessage"]      = WarningMessage;
      name_map["DebugMessage"]        = DebugMessage;
      name_map["VerboseDebugMessage"] = VerboseDebugMessage;
      name_map["EveryMessage"]        = EveryMessage;
      name_map["*"]                   = EveryMessage;
    }
    MapType::const_iterator level = name_map.find(name);
    if (level != name_map.end())
      return level->second;

    return boost::lexical_cast<int32>(name);
  }
}

void vw::parse_config_file(const char* fn, vw::Settings& settings) {
  ifstream file(fn);

  if (!file.is_open())
    vw_throw(IOErr() <<  "Could not open logfile: " << fn);

  parse_config(file, settings);
}

void vw::parse_config(basic_istream<char>& stream, vw::Settings& settings) {

  // DO NOT try to log with vw_log! It will cause a deadlock because we're
  // holding locks inside reload_config, and the loggers will call
  // reload_config first.
  po::parsed_options opts(0);

  try {
     opts = po::parse_config_file( stream, po::options_description(), true );
  } catch (const po::invalid_syntax& e) {
    cerr << "Could not parse config file. Ignoring. (" << e.msg << " near \"" << e.tokens << "\")" << endl;
  }

  shared_ptr<LogInstance> current_log;
  string current_logname = "console";

  BOOST_FOREACH( const po::option& o, opts.options ) {

    try {
      if (o.string_key == "general.default_num_threads")
        settings.set_default_num_threads(boost::lexical_cast<int>(o.value[0]));
      else if (o.string_key == "general.system_cache_size")
        settings.set_system_cache_size(boost::lexical_cast<size_t>(o.value[0]));
      else if (o.string_key.compare(0, 8, "logfile ") == 0) {
        size_t sep = o.string_key.find_last_of('.');
        assert(sep != string::npos);

        string logname = o.string_key.substr(8, sep-8);
        string level_s = o.string_key.substr(sep+1);
        string domain  = o.value[0];

        cerr << "Logname[" << logname << "] level_s[" << level_s << "] domain[" << domain << "]" << endl;

        if (logname.empty() || level_s.empty() || domain.empty())
          continue;

        if (logname != current_logname) {
          current_logname = logname;
          if (current_logname == "console")
            current_log.reset();
          else {
            current_log = shared_ptr<LogInstance>( new LogInstance(logname) );
            vw_log().add(current_log);
          }
        }

        int32 level = name2level(level_s);

        cerr << "Adding rule: " << level << " " << domain << "\n";
        if (current_log)
          current_log->rule_set().add_rule(level, domain);
        else
          vw_log().console_log().rule_set().add_rule(level, domain);
      }
      else {
        continue;
      }
    } catch (const boost::bad_lexical_cast& e) {
      cerr << "Bad line in config file near " << o.string_key << ". skipping." << endl;
    }
  }

}
