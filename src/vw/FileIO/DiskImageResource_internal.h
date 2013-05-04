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


/// \file FileIO/DiskImageResource_internal.h
///
/// A header for internal use only that allows access to the extension list
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL_H__
#define __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL_H__

#include <string>
#include <set>
#include <boost/function.hpp>

namespace vw {
namespace internal {
  typedef boost::function<void (std::string const&)> ExtTestFunction;
  void foreach_ext(std::string const& prefix, ExtTestFunction const& callback,
                  std::set<std::string> const& exclude = std::set<std::string>() );
}} // namespace vw::internal

#endif // __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL__
