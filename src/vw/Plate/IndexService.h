// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_INDEX_SERVICE_H__
#define __VW_PLATE_INDEX_SERVICE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/Protobuffers.pb.h>

#include <google/protobuf/service.h>

namespace vw {
namespace platefile {

  class IndexServiceImpl : public IndexService {
    std::string m_root_directory;
    Index m_index;


    // Private methods
    std::vector<std::string> plate_filenames(std::string const& root_directory);

  public:

    IndexServiceImpl(std::string root_directory);

    void OpenRequest(::google::protobuf::RpcController* controller,
                     const ::vw::platefile::IndexOpenRequest* request,
                     ::vw::platefile::IndexOpenReply* response,
                     ::google::protobuf::Closure* done);

    void ReadRequest(::google::protobuf::RpcController* controller,
                     const ::vw::platefile::IndexReadRequest* request,
                     ::vw::platefile::IndexReadReply* response,
                     ::google::protobuf::Closure* done);      
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_INDEX_SERVICE_H__
