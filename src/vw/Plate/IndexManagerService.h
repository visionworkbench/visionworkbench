// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_INDEX_MANAGER_SERVICE_H__
#define __VW_PLATE_INDEX_MANAGER_SERVICE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <google/protobuf/service.h>

namespace vw {
namespace platefile {

  class IndexManagerServiceImpl : public IndexManagerService {

  public:

    IndexManagerServiceImpl();

    virtual void ListRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexListRequest* request,
                             ::vw::platefile::IndexListReply* response,
                             ::google::protobuf::Closure* done);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_INDEX_SERVICE_H__
