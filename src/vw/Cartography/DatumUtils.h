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

#ifndef __VW_CARTOGRAPHY_DATUM_UTILS_H__
#define __VW_CARTOGRAPHY_DATUM_UTILS_H__

#include <vw/Cartography/Datum.h>

/// \file DatumUtils.h Tools for working with Datum objects

namespace vw {

// This can catch user mistakes
void checkDatumConsistency(vw::cartography::Datum const& datum1,
                           vw::cartography::Datum const& datum2,
                           bool warn_only);

} // namespace vw

#endif // __VW_CARTOGRAPHY_DATUM_UTILS_H__
