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


#include <vw/Plate/Exception.h>

void vw::platefile::throw_rpc_error(const RpcErrorMsg& e) {
  switch (e.code()) {
    case RpcErrorMsg::SUCCESS:
      return;
    case RpcErrorMsg::UNKNOWN:
      vw_throw(RpcErr() << "Unrecognized Error: " << e.msg());
    case RpcErrorMsg::REMOTE_ERROR:
      vw_throw(RpcErr() << "Server Error: " << e.msg());
    case RpcErrorMsg::PLATE_UNKNOWN:
      vw_throw(PlatefileErr() << e.msg());
    case RpcErrorMsg::PLATE_TILE_NOT_FOUND:
      vw_throw(TileNotFoundErr() << e.msg());
    case RpcErrorMsg::PLATE_INVALID:
      vw_throw(InvalidPlatefileErr() << e.msg());
    case RpcErrorMsg::PLATE_CREATION_FAILED:
      vw_throw(PlatefileCreationErr() << e.msg());
    case RpcErrorMsg::PLATE_BLOB_IO:
      vw_throw(BlobIoErr() << e.msg());
    case RpcErrorMsg::PLATE_BLOB_LIMIT:
      vw_throw(BlobLimitErr() << e.msg());
  }
}
