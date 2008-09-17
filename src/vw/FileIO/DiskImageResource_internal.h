// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file FileIO/DiskImageResource_internal.h
///
/// A header for internal use only that allows access to the extension list
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL_H__
#define __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL_H__

#include <string>
#include <set>

namespace vw {
namespace internal {

  typedef void (*ExtTestFunction)(std::string const& fn);
  void foreach_ext(std::string const& fn, ExtTestFunction func, std::set<std::string> const& exclude = std::set<std::string>() );

}} // namespace vw::internal

#endif // __VW_FILEIO_DISKIMAGERESOURCE_INTERNAL__
