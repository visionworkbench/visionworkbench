// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/IndexService.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/IndexPage.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Log.h>

using namespace vw::platefile;

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------------
//                                 PRIVATE METHODS
// ----------------------------------------------------------------------------------

std::vector<std::string> IndexServiceImpl::glob_plate_filenames(std::string const& root_directory)  {

  std::vector<std::string> result;

  if ( !fs::exists( root_directory ) )
    vw_throw(IOErr() << "Could not access the root directory: \"" << root_directory << "\"");

  // Create a regular expression for matching the pattern 'plate_<number>.blob"
  boost::regex re;
  re.assign(".*\\.plate", boost::regex_constants::icase);

  fs::directory_iterator end_itr; // default construction yields past-the-end

  // Iterate through the files in the platefile directory and return
  // any that match the regex above.
  for ( fs::directory_iterator itr( root_directory ); itr != end_itr; ++itr ) {
    if (boost::regex_match(itr->leaf(), re))
      result.push_back(itr->leaf());
  }

  return result;
}

IndexServiceImpl::IndexServiceRecord IndexServiceImpl::add_index(std::string root_directory,
                                                                 std::string plate_filename,
                                                                 boost::shared_ptr<Index> index) {

    // Build up an IndexServiceRecord
    IndexServiceRecord rec;
    rec.short_plate_filename = plate_filename;
    rec.full_plate_filename = root_directory + "/" + plate_filename;
    rec.index = index;

    // Store the record in a std::map by platefile_id
    m_indices[index->index_header().platefile_id()]  = rec;
    return rec;
}

/// Fetch an IndexServiceRecord for a given platefile_id, or throw an
/// exception if no record is found.  This is the first step in most
/// of the RPC handlers below, so the code is factored out here for
/// convenience.
IndexServiceImpl::IndexServiceRecord IndexServiceImpl::get_index_record_for_platefile_id(int platefile_id) {

  if (m_indices.find(platefile_id) == m_indices.end())
    vw_throw(InvalidPlatefileErr() << "No platefile matching this platefile_id could be found.");

  return m_indices[platefile_id];

}

// ----------------------------------------------------------------------------------
//                                 PUBLIC METHODS
// ----------------------------------------------------------------------------------

IndexServiceImpl::IndexServiceImpl(std::string root_directory) :
  m_root_directory( fs::system_complete(root_directory).string() ) {

  // Search for all platefiles in the given root_directory.  A
  // platefile is any directory ending in *.plate.
  std::vector<std::string> platefiles = glob_plate_filenames(m_root_directory);
  if (platefiles.size() < 1)
    vw_out(InfoMessage, "plate:index_service") << "Warning: could not find any platefiles in the root directory.";

  for (unsigned i = 0 ; i < platefiles.size(); ++i) {
    // Open each new platefile.  This will return a LocalIndex which we store using add_index()
    boost::shared_ptr<Index> idx = Index::construct_open(Url(m_root_directory + "/" + platefiles[i]));
    this->add_index(m_root_directory, platefiles[i], idx);
  }
}

void IndexServiceImpl::sync() {
  for (index_list_type::iterator iter = m_indices.begin(); iter != m_indices.end(); ++iter) {
    vw_out() << "\t--> Syncing index for " << iter->second.short_plate_filename
             << " to disk.\n";
    iter->second.index->sync();
  }
}

void IndexServiceImpl::OpenRequest(::google::protobuf::RpcController* /*controller*/,
                                   const IndexOpenRequest* request,
                                   IndexOpenReply* response,
                                   ::google::protobuf::Closure* done) {

  index_list_type::iterator iter = m_indices.begin();
  while (iter != m_indices.end()) {
    if ( (*iter).second.short_plate_filename == request->plate_name() ) {

      // If we find a matching short_plate_filename in m_indices, then
      // we return the index header.
      response->set_short_plate_filename( (*iter).second.short_plate_filename );
      response->set_full_plate_filename( (*iter).second.full_plate_filename );
      *(response->mutable_index_header()) = (*iter).second.index->index_header();
      done->Run();
      return;

    }
    ++iter;
  }

  // If we reach this point, then there is no matching plate file
  // being tracked by this IndexService.  Return an error.
  vw_throw(InvalidPlatefileErr() << "No platefile matching this plate_name could be found.");
}

void IndexServiceImpl::CreateRequest(::google::protobuf::RpcController* /*controller*/,
                                     const IndexCreateRequest* request,
                                     IndexOpenReply* response,
                                     ::google::protobuf::Closure* done) {

  // Try OPEN: the index may already exist and be open.  If so, we do
  // nothing but return the necessary info in the IndexOpenReply.
  index_list_type::iterator iter = m_indices.begin();
  while (iter != m_indices.end()) {
    if ( (*iter).second.short_plate_filename == request->plate_name() ) {

      // If we find a matching short_plate_filename in m_indices, then
      // we return the index header.
      response->set_short_plate_filename( (*iter).second.short_plate_filename );
      response->set_full_plate_filename( (*iter).second.full_plate_filename );
      *(response->mutable_index_header()) = (*iter).second.index->index_header();
      done->Run();
      return;

    }
    ++iter;
  }


  // If the index is not already open...
  std::string url = m_root_directory + "/" + request->plate_name();
  boost::shared_ptr<Index> idx;
  try {
    if (!exists( fs::path( url, fs::native ) ) ) {

      // CREATE: If there was no index already open by that name, and
      // no platefile on disk, then we try to create it.
      fs::create_directory(url);
      idx = Index::construct_create(url, request->index_header());

    } else {

      // REOPEN: If a platefile actually does appear to exist on disk,
      // but we had not previously opened it, we open it here.
      idx = Index::construct_open(url);

    }
  } catch (IOErr &e) {
    vw_throw(PlatefileCreationErr() << e.what());
    return; // never reached
  }

  IndexServiceRecord rec = this->add_index(m_root_directory, request->plate_name(), idx);
  response->set_short_plate_filename( rec.short_plate_filename );
  response->set_full_plate_filename( rec.full_plate_filename );
  *(response->mutable_index_header()) = rec.index->index_header();
  done->Run();
}

void IndexServiceImpl::InfoRequest(::google::protobuf::RpcController* /*controller*/,
                                   const IndexInfoRequest* request,
                                   IndexInfoReply* response,
                                   ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  response->set_short_plate_filename(rec.short_plate_filename);
  response->set_full_plate_filename(rec.full_plate_filename);
  *(response->mutable_index_header()) = rec.index->index_header();
  done->Run();
}

void IndexServiceImpl::ListRequest(::google::protobuf::RpcController* /*controller*/,
                                   const IndexListRequest* request,
                                   IndexListReply* response,
                                   ::google::protobuf::Closure* done) {

  //  std::cout << "RECEIVED ListRequest message : " << request->DebugString() << "\n";

  for (index_list_type::const_iterator i = m_indices.begin(), end = m_indices.end(); i != end; ++i) {

    IndexServiceRecord rec = i->second;

#define should_filter_out(field) \
    request->has_ ## field() && rec.index->index_header().has_ ## field() && request->field() != rec.index->index_header().field()

    if (should_filter_out(tile_filetype))
      continue;

    if (should_filter_out(pixel_format))
      continue;

    if (should_filter_out(channel_type))
      continue;

    if (should_filter_out(type))
      continue;

#undef should_filter_out

    response->add_platefile_names(rec.short_plate_filename);
  }

  //  std::cout << "SENDING ListRequest response : " << response->DebugString() << "\n";

  done->Run();
}

void IndexServiceImpl::PageRequest(::google::protobuf::RpcController* /*controller*/,
                                   const IndexPageRequest* request,
                                   IndexPageReply* response,
                                   ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  boost::shared_ptr<IndexPage> page = rec.index->page_request(request->col(),
                                                              request->row(),
                                                              request->level());

  std::ostringstream ostr;
  page->serialize(ostr);
  response->set_page_bytes(ostr.str().c_str(), ostr.str().size());
  done->Run();
}

void IndexServiceImpl::ReadRequest(::google::protobuf::RpcController* /*controller*/,
                                   const IndexReadRequest* request,
                                   IndexReadReply* response,
                                   ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  IndexRecord record = rec.index->read_request(request->col(),
                                               request->row(),
                                               request->level(),
                                               request->transaction_id(),
                                               request->exact_transaction_match());
  *(response->mutable_index_record()) = record;
  done->Run();
}

void IndexServiceImpl::WriteRequest(::google::protobuf::RpcController* /*controller*/,
                                    const IndexWriteRequest* request,
                                    IndexWriteReply* response,
                                    ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  uint64 size;
  int blob_id = rec.index->write_request(size);
  response->set_blob_id(blob_id);
  response->set_size(boost::numeric_cast<int32>(size));
  done->Run();
}

void IndexServiceImpl::WriteUpdate(::google::protobuf::RpcController* /*controller*/,
                                   const ::vw::platefile::IndexWriteUpdate* request,
                                   ::vw::platefile::RpcNullMsg* /*response*/,
                                   ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  rec.index->write_update(request->header(), request->record());

  // This message has no response
  done->Run();
}

void IndexServiceImpl::MultiWriteUpdate(::google::protobuf::RpcController* /*controller*/,
                                        const ::vw::platefile::IndexMultiWriteUpdate* request,
                                        ::vw::platefile::RpcNullMsg* /*response*/,
                                        ::google::protobuf::Closure* done) {

  // MultiWrite updates are packetized.  We iterate over them here.
  for (int i = 0; i < request->write_updates().size(); ++i) {
    IndexWriteUpdate update = request->write_updates().Get(i);

    // Fetch the index service record
    IndexServiceRecord rec = get_index_record_for_platefile_id(update.platefile_id());

    // Access the data in the index.  Return the data on success, or
    // notify the remote client of our failure if we did not succeed.
    rec.index->write_update(update.header(), update.record());
  }

  // This message has no response
  done->Run();
}

void IndexServiceImpl::WriteComplete(::google::protobuf::RpcController* /*controller*/,
                                     const IndexWriteComplete* request,
                                     RpcNullMsg* /*response*/,
                                     ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  rec.index->write_complete(request->blob_id(), request->blob_offset());

  // This message has no response
  done->Run();
}


void IndexServiceImpl::TransactionRequest(::google::protobuf::RpcController* /*controller*/,
                                          const IndexTransactionRequest* request,
                                          IndexTransactionReply* response,
                                          ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  int transaction_id = rec.index->transaction_request(request->description(),
                                                      request->transaction_id_override());
  response->set_transaction_id(transaction_id);

  // std::cout << "\n\t--> [ Platefile " << request->platefile_id() << " ] : Transaction "
  //           << transaction_id << " started.\n";

  done->Run();
}

void IndexServiceImpl::TransactionComplete(::google::protobuf::RpcController* /*controller*/,
                                           const IndexTransactionComplete* request,
                                           RpcNullMsg* /*response*/,
                                           ::google::protobuf::Closure* done) {

  // std::cout << "\n\t--> [ Platefile " << request->platefile_id() << " ] : Transaction "
  //           << request->transaction_id() << " complete.\n";

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  rec.index->transaction_complete(request->transaction_id(),
                                  request->update_read_cursor());

  // This message has no response.
  done->Run();
}

void IndexServiceImpl::TransactionFailed(::google::protobuf::RpcController* /*controller*/,
                                         const IndexTransactionFailed* request,
                                         RpcNullMsg* /*response*/,
                                         ::google::protobuf::Closure* done) {

  // std::cout << "\n\t--> [ Platefile " << request->platefile_id() << " ] : Transaction "
  //           << request->transaction_id() << " failed.\n";

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  rec.index->transaction_failed(request->transaction_id());

  // This message has no response.
  done->Run();
}

void IndexServiceImpl::TransactionCursor(::google::protobuf::RpcController* /*controller*/,
                                         const IndexTransactionCursorRequest* request,
                                         IndexTransactionCursorReply* response,
                                         ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  int transaction_id = rec.index->transaction_cursor();
  response->set_transaction_id(transaction_id);
  done->Run();
}

void IndexServiceImpl::NumLevelsRequest(::google::protobuf::RpcController* /*controller*/,
                                        const IndexNumLevelsRequest* request,
                                        IndexNumLevelsReply* response,
                                        ::google::protobuf::Closure* done) {

  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  response->set_num_levels(rec.index->num_levels());
  done->Run();
}

void IndexServiceImpl::LogRequest(::google::protobuf::RpcController* /*controller*/,
                                  const ::vw::platefile::IndexLogRequest* request,
                                  ::vw::platefile::RpcNullMsg* /*response*/,
                                  ::google::protobuf::Closure* done) {
  // Fetch the index service record
  IndexServiceRecord rec = get_index_record_for_platefile_id(request->platefile_id());
  rec.index->log(request->message());

  // This message has no response
  done->Run();
}

void IndexServiceImpl::TestRequest(::google::protobuf::RpcController* /*controller*/,
                                    const IndexTestRequest* request,
                                    IndexTestReply* response,
                                    ::google::protobuf::Closure* done) {

  // access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  response->set_value(request->value());
  done->Run();
}

