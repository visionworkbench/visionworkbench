// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_INDEX_SERVICE_H__
#define __VW_PLATE_INDEX_SERVICE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <google/protobuf/service.h>

namespace vw {
namespace platefile {

  class IndexServiceImpl : public IndexService {


    struct IndexServiceRecord {
      std::string short_plate_filename;
      std::string full_plate_filename;
      boost::shared_ptr<Index> index;
    };

    std::string m_root_directory;

    typedef std::map<int32, IndexServiceRecord> index_list_type;
    index_list_type m_indices;

    // Private methods
    std::vector<std::string> glob_plate_filenames(std::string const& root_directory);
    IndexServiceRecord add_index(std::string root_directory, std::string plate_filename,
                                 boost::shared_ptr<Index> index);
      
    /// Fetch an IndexServiceRecord for a given platefile_id, or throw an
    /// exception if no record is found.
    IndexServiceRecord get_index_record_for_platefile_id(int platefile_id);

  public:

    IndexServiceImpl(std::string root_directory);

    virtual void OpenRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexOpenRequest* request,
                             ::vw::platefile::IndexOpenReply* response,
                             ::google::protobuf::Closure* done);

    virtual void CreateRequest(::google::protobuf::RpcController* controller,
                               const ::vw::platefile::IndexCreateRequest* request,
                               ::vw::platefile::IndexOpenReply* response,
                               ::google::protobuf::Closure* done);

    virtual void InfoRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexInfoRequest* request,
                             ::vw::platefile::IndexInfoReply* response,
                             ::google::protobuf::Closure* done);

    virtual void ListRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexListRequest* request,
                             ::vw::platefile::IndexListReply* response,
                             ::google::protobuf::Closure* done);

    virtual void ReadRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexReadRequest* request,
                             ::vw::platefile::IndexReadReply* response,
                             ::google::protobuf::Closure* done);      

    virtual void MultiReadRequest(::google::protobuf::RpcController* controller,
                                  const ::vw::platefile::IndexMultiReadRequest* request,
                                  ::vw::platefile::IndexMultiReadReply* response,
                                  ::google::protobuf::Closure* done);      

    virtual void WriteRequest(::google::protobuf::RpcController* controller,
                              const ::vw::platefile::IndexWriteRequest* request,
                              ::vw::platefile::IndexWriteReply* response,
                              ::google::protobuf::Closure* done);      

    virtual void WriteComplete(::google::protobuf::RpcController* controller,
                               const ::vw::platefile::IndexWriteComplete* request,
                               ::vw::platefile::RpcNullMessage* response,
                               ::google::protobuf::Closure* done);      

    virtual void TransactionRequest(::google::protobuf::RpcController* controller,
                                    const ::vw::platefile::IndexTransactionRequest* request,
                                    ::vw::platefile::IndexTransactionReply* response,
                                    ::google::protobuf::Closure* done);      

    virtual void TransactionComplete(::google::protobuf::RpcController* controller,
                                     const ::vw::platefile::IndexTransactionComplete* request,
                                     ::vw::platefile::RpcNullMessage* response,
                                     ::google::protobuf::Closure* done);      

    virtual void TransactionFailed(::google::protobuf::RpcController* controller,
                                   const ::vw::platefile::IndexTransactionFailed* request,
                                   ::vw::platefile::RpcNullMessage* response,
                                   ::google::protobuf::Closure* done);      

    virtual void TransactionCursor(::google::protobuf::RpcController* controller,
                                   const ::vw::platefile::IndexTransactionCursorRequest* request,
                                   ::vw::platefile::IndexTransactionCursorReply* response,
                                   ::google::protobuf::Closure* done);      

    virtual void ValidTiles(::google::protobuf::RpcController* controller,
                            const ::vw::platefile::IndexValidTilesRequest* request,
                            ::vw::platefile::IndexValidTilesReply* response,
                            ::google::protobuf::Closure* done);      
    
    virtual void NumLevelsRequest(::google::protobuf::RpcController* controller,
                                  const ::vw::platefile::IndexNumLevelsRequest* request,
                                  ::vw::platefile::IndexNumLevelsReply* response,
                                  ::google::protobuf::Closure* done);      

    // A simple message that echos back the value that was sent.
    virtual void TestRequest(::google::protobuf::RpcController* controller,
                             const ::vw::platefile::IndexTestRequest* request,
                             ::vw::platefile::IndexTestReply* response,
                             ::google::protobuf::Closure* done);      
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_INDEX_SERVICE_H__
