// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
