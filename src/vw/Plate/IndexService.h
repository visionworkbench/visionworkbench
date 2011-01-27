// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_INDEX_SERVICE_H__
#define __VW_PLATE_INDEX_SERVICE_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexService.pb.h>
#include <boost/shared_ptr.hpp>

namespace vw {
namespace platefile {

  class Index;

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
    IndexServiceRecord* add_index(std::string root_directory, std::string plate_filename,
                                  boost::shared_ptr<Index> index);

    IndexServiceRecord* find_id(int32 platefile_id);
    IndexServiceRecord  find_id_throw(int32 platefile_id);
    IndexServiceRecord* find_name(const std::string& name);

  public:

    IndexServiceImpl(std::string root_directory);

    void sync();

    virtual void OpenRequest(::google::protobuf::RpcController* controller,
                             const IndexOpenRequest* request,
                             IndexOpenReply* response,
                             ::google::protobuf::Closure* done);

    virtual void CreateRequest(::google::protobuf::RpcController* controller,
                               const IndexCreateRequest* request,
                               IndexOpenReply* response,
                               ::google::protobuf::Closure* done);

    virtual void InfoRequest(::google::protobuf::RpcController* controller,
                             const IndexInfoRequest* request,
                             IndexInfoReply* response,
                             ::google::protobuf::Closure* done);

    virtual void ListRequest(::google::protobuf::RpcController* controller,
                             const IndexListRequest* request,
                             IndexListReply* response,
                             ::google::protobuf::Closure* done);

    virtual void PageRequest(::google::protobuf::RpcController* controller,
                             const IndexPageRequest* request,
                             IndexPageReply* response,
                             ::google::protobuf::Closure* done);

    virtual void ReadRequest(::google::protobuf::RpcController* controller,
                             const IndexReadRequest* request,
                             IndexReadReply* response,
                             ::google::protobuf::Closure* done);

    virtual void WriteRequest(::google::protobuf::RpcController* controller,
                              const IndexWriteRequest* request,
                              IndexWriteReply* response,
                              ::google::protobuf::Closure* done);

    virtual void WriteUpdate(::google::protobuf::RpcController* controller,
                             const IndexWriteUpdate* request,
                             RpcNullMsg* response,
                             ::google::protobuf::Closure* done);

    // Like WriteUpdate, but packetizes updates.
    virtual void MultiWriteUpdate(::google::protobuf::RpcController* controller,
                             const IndexMultiWriteUpdate* request,
                             RpcNullMsg* response,
                             ::google::protobuf::Closure* done);

    virtual void WriteComplete(::google::protobuf::RpcController* controller,
                               const IndexWriteComplete* request,
                               RpcNullMsg* response,
                               ::google::protobuf::Closure* done);

    virtual void TransactionRequest(::google::protobuf::RpcController* controller,
                                    const IndexTransactionRequest* request,
                                    IndexTransactionReply* response,
                                    ::google::protobuf::Closure* done);

    virtual void TransactionComplete(::google::protobuf::RpcController* controller,
                                     const IndexTransactionComplete* request,
                                     RpcNullMsg* response,
                                     ::google::protobuf::Closure* done);

    virtual void TransactionFailed(::google::protobuf::RpcController* controller,
                                   const IndexTransactionFailed* request,
                                   RpcNullMsg* response,
                                   ::google::protobuf::Closure* done);

    virtual void TransactionCursor(::google::protobuf::RpcController* controller,
                                   const IndexTransactionCursorRequest* request,
                                   IndexTransactionCursorReply* response,
                                   ::google::protobuf::Closure* done);

    virtual void NumLevelsRequest(::google::protobuf::RpcController* controller,
                                  const IndexNumLevelsRequest* request,
                                  IndexNumLevelsReply* response,
                                  ::google::protobuf::Closure* done);

    virtual void LogRequest(::google::protobuf::RpcController* controller,
                            const IndexLogRequest* request,
                            RpcNullMsg* response,
                            ::google::protobuf::Closure* done);

    // A simple message that echos back the value that was sent.
    virtual void TestRequest(::google::protobuf::RpcController* controller,
                             const IndexTestRequest* request,
                             IndexTestReply* response,
                             ::google::protobuf::Closure* done);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_INDEX_SERVICE_H__
