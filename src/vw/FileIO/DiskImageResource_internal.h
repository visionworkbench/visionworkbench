// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
