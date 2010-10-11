// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Shaders.h

#ifndef VW_GPU_SHADERS_H
#define VW_GPU_SHADERS_H

#include <map>
#include <string>

namespace vw { namespace GPU {

extern std::map<std::string, const char*> standard_shaders_map;

void init_standard_shaders();

#endif

 } }
