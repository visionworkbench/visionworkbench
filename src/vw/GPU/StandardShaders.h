
// StandardShaders.h

#ifndef StandardShaders_H
#define StandardShaders_H

#include <map>
#include <string>

namespace vw { namespace GPU {

extern std::map<std::string, char*> standard_shaders_map;

void init_standard_shaders();

#endif

 } } 
