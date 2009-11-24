// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/IndexService.h>
using namespace vw::platefile;

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

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
    rec.index_header = index->index_header();
    rec.index = index;
    
    // Store the record in a std::map by platefile_id
    m_indices[index->index_header().platefile_id()]  = rec;
    return rec;
}

IndexServiceImpl::IndexServiceImpl(std::string root_directory) : 
  m_root_directory(root_directory) {
  
  std::cout << "Starting Index Service\n";

  // Search for all platefiles in the given root_directory.  A
  // platefile is any directory ending in *.plate.
  std::vector<std::string> platefiles = glob_plate_filenames(fs::system_complete(root_directory));
  if (platefiles.size() < 1) 
    vw_throw(ArgumentErr() << "Error: could not find any platefiles in the root directory.");
  
  for (unsigned i = 0 ; i < platefiles.size(); ++i) {

    // Open each new platefile.  This will return a LocalIndex which we store using add_index()
    boost::shared_ptr<Index> idx = Index::construct_open(root_directory + "/" + platefiles[i]);
    this->add_index(root_directory, platefiles[i], idx);
  }    

}

void IndexServiceImpl::OpenRequest(::google::protobuf::RpcController* controller,
                                   const IndexOpenRequest* request,
                                   IndexOpenReply* response,
                                   ::google::protobuf::Closure* done) {

  index_list_type::iterator iter = m_indices.begin();
  while (iter != m_indices.end()) {
    if ( (*iter).second.short_plate_filename == request->plate_name() ) {

      // If we find a matching short_plate_filename in m_indices, then
      // we return the index header.
      *(response->mutable_index_header()) = (*iter).second.index_header;
      done->Run();    
      return;

    }
    ++iter;
  }

  // If we reach this point, then there is no matching plate file
  // being tracked by this IndexService.  Return an error.
  
  controller->SetFailed("No platefile matching this plate_name could be found.");
  done->Run();    
}

void IndexServiceImpl::CreateRequest(::google::protobuf::RpcController* controller,
                                     const IndexCreateRequest* request,
                                     IndexOpenReply* response,
                                     ::google::protobuf::Closure* done) {

  std::string url = m_root_directory + "/" + request->plate_name();
  if( exists( fs::path( url, fs::native ) ) ) {
    controller->SetFailed("There is already an existing platefile with that name.");
    done->Run();
    return;
  }

  boost::shared_ptr<Index> idx;
  try {
    fs::create_directory(url);
    idx = Index::construct_create(url, request->index_header());
  } catch (IOErr &e) {
    controller->SetFailed(e.what());
    done->Run();
    return;
  }

  IndexServiceRecord rec = this->add_index(m_root_directory, request->plate_name(), idx);
  *(response->mutable_index_header()) = rec.index_header;
  done->Run();
}

void IndexServiceImpl::InfoRequest(::google::protobuf::RpcController* controller,
                                   const IndexInfoRequest* request,
                                   IndexInfoReply* response,
                                   ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }

  IndexServiceRecord rec = m_indices[request->platefile_id()];
  response->set_short_plate_filename(rec.short_plate_filename);
  response->set_full_plate_filename(rec.full_plate_filename);
  *(response->mutable_index_header()) = rec.index_header;
  done->Run();
}

void IndexServiceImpl::ReadRequest(::google::protobuf::RpcController* controller,
                                   const IndexReadRequest* request,
                                   IndexReadReply* response,
                                   ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    IndexRecord record = rec.index->read_request(request->col(), 
                                                 request->row(), 
                                                 request->depth(), 
                                                 request->transaction_id());
    *(response->mutable_index_record()) = record;
  } catch (TileNotFoundErr &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}


void IndexServiceImpl::WriteRequest(::google::protobuf::RpcController* controller,
                                    const IndexWriteRequest* request,
                                    IndexWriteReply* response,
                                    ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    int blob_id = rec.index->write_request(request->size());
    response->set_blob_id(blob_id);
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}

void IndexServiceImpl::WriteComplete(::google::protobuf::RpcController* controller,
                                     const IndexWriteComplete* request,
                                     RpcNullMessage* response,
                                     ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    rec.index->write_complete(request->header(), request->record());
    // This message has no response
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}

void IndexServiceImpl::TransactionRequest(::google::protobuf::RpcController* controller,
                                          const IndexTransactionRequest* request,
                                          IndexTransactionReply* response,
                                          ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    int transaction_id = rec.index->transaction_request(request->description());
    response->set_transaction_id(transaction_id);
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}   

void IndexServiceImpl::TransactionComplete(::google::protobuf::RpcController* controller,
                                           const IndexTransactionComplete* request,
                                           RpcNullMessage* response,
                                           ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    rec.index->transaction_complete(request->transaction_id());
    // This message has no response.
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}

void IndexServiceImpl::TransactionCursor(::google::protobuf::RpcController* controller,
                                         const IndexTransactionCursorRequest* request,
                                         IndexTransactionCursorReply* response,
                                         ::google::protobuf::Closure* done) {

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];
    
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    int transaction_id = rec.index->transaction_cursor();
    response->set_transaction_id(transaction_id);
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}

void IndexServiceImpl::DepthRequest(::google::protobuf::RpcController* controller,
                                    const IndexDepthRequest* request,
                                    IndexDepthReply* response,
                                    ::google::protobuf::Closure* done) {
  std::cout << "Processing depth request....\n";

  if (m_indices.find(request->platefile_id()) == m_indices.end()) {
    controller->SetFailed("Invalid Platefile ID.");
    done->Run();
    return;
  }
  IndexServiceRecord rec = m_indices[request->platefile_id()];

  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    response->set_depth(rec.index->max_depth());
  } catch (vw::Exception &e) {
    controller->SetFailed(e.what());
  }
  done->Run();
}

