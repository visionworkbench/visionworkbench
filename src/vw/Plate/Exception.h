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


// \file Exception.h
//
// Exceptions used in the platefile system.
//
#ifndef __VW_PLATE_EXCEPTION_H__
#define __VW_PLATE_EXCEPTION_H__

#include <vw/Core/Exception.h>
#include <vw/Plate/Rpc.pb.h>

namespace vw {
namespace platefile {

  // Three types of exceptions we care about:
  //   1) In-band client errors (eg, asked for an invalid platefile)
  //   2) In-band server errors (eg, server ran out of memory)
  //   3) transport errors (eg, channel closed)

  // The in-band errors are represented as RpcErr.
  //  The in-band client errors are represented as PlatefileErr
  // The out-of-band errors are represented as subclasses of NetworkErr

  VW_DEFINE_EXCEPTION(RpcErr, IOErr);
  VW_DEFINE_EXCEPTION(NetworkErr, IOErr);

  VW_DEFINE_EXCEPTION_EXT(PlatefileErr, RpcErr) {
    VW_EXCEPTION_API(PlatefileErr);
    virtual RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_UNKNOWN; }
  };

  // This exception is thrown by the Tree and Index classes whenever
  // a tile is requested that does not exist.  It is frequently
  // caught by higher level classes like PlateFile when they are
  // trying to determine whether a tile exists or not.
  VW_DEFINE_EXCEPTION_EXT(TileNotFoundErr, PlatefileErr) {
    VW_EXCEPTION_API(TileNotFoundErr);
    RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_TILE_NOT_FOUND; }
  };

  // This exception is thrown by the IndexService when there is an
  // attempt to access or operate on a platefile that has not been
  // opened, and is not being tracked by the system.
  VW_DEFINE_EXCEPTION_EXT(InvalidPlatefileErr, PlatefileErr) {
    VW_EXCEPTION_API(InvalidPlatefileErr);
    RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_INVALID; }
  };

  // This exception is thrown whenever the maximum number of blobs is
  // reached.  This shouldn't ever really happen in practice unless
  // something is wrong.
  VW_DEFINE_EXCEPTION_EXT(BlobLimitErr, PlatefileErr) {
    VW_EXCEPTION_API(BlobLimitErr);
    RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_BLOB_LIMIT; }
  };

  // This exception is thrown when something bad happens while reading or
  // writing a blob.
  VW_DEFINE_EXCEPTION_EXT(BlobIoErr, PlatefileErr) {
    VW_EXCEPTION_API(BlobIoErr);
    RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_BLOB_IO; }
  };

  // This exception is thrown by the IndexService when an error
  // occurs while attempting to create a new platefile.
  VW_DEFINE_EXCEPTION_EXT(PlatefileCreationErr, PlatefileErr) {
    VW_EXCEPTION_API(PlatefileCreationErr);
    RpcErrorMsg::Code code() const { return RpcErrorMsg::PLATE_CREATION_FAILED; }
  };

  void throw_rpc_error(const RpcErrorMsg& e);

}} // namespace vw::platefile

#endif // __VW_PLATE_TREE_H__
