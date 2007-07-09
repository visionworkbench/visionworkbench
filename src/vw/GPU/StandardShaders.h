
// StandardShaders.h

#ifndef StandardShaders_H
#define StandardShaders_H

#include <map>
#include <string>

namespace vw { namespace GPU {

extern std::map<std::string, char*> standard_shader_map;

void init_standard_shaders();

#endif

 } } 
