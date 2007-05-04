%module core

%{
#include <vw/Core.h>
%}

namespace vw {

  enum MessageLevel {
    ErrorMessage = 0,
    WarningMessage = 10,
    InfoMessage = 20,
    DebugMessage = 30,
    VerboseDebugMessage = 40
  };

  void set_debug_level( int level );

} // namespace vw
