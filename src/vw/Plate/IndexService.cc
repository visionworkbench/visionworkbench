// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/IndexService.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


vw::platefile::IndexServiceImpl::IndexServiceImpl(std::string root_directory) : 
  m_root_directory(root_directory),
  m_index("/opt/local/apache2/htdocs/platefiles/TrueMarble.16km.2700x1350_toast.plate") {

  std::cout << "Starting Index Service\n";
  std::vector<std::string> platefiles = plate_filenames(root_directory);
  if (platefiles.size() < 1) {
    std::cout << "Error: could not find any platefiles in the root directory."
              << "\nExiting.\n\n";
    exit(1);
  }
  
  for (int i = 0 ; i<platefiles.size(); ++i) {
    std::cout << "\t--> Found platefile: " 
              << root_directory << "/" << platefiles[i] << "\n";
  }    

  // TODO: finish this function
  // Index idx(root_directory + "/" + platefiles[0]);
}



std::vector<std::string> vw::platefile::IndexServiceImpl::plate_filenames(std::string const& root_directory)  {

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


void vw::platefile::IndexServiceImpl::OpenRequest(::google::protobuf::RpcController* controller,
                                   const ::vw::platefile::IndexOpenRequest* request,
                                   ::vw::platefile::IndexOpenReply* response,
                                   ::google::protobuf::Closure* done) {

  
  // These are placeholders for now
  response->set_platefile_id(23);
  response->set_secret(11);
  
  done->Run();
}

void vw::platefile::IndexServiceImpl::ReadRequest(::google::protobuf::RpcController* controller,
                                   const ::vw::platefile::IndexReadRequest* request,
                                   ::vw::platefile::IndexReadReply* response,
                                   ::google::protobuf::Closure* done) {
  
  vw_out(InfoMessage, "platefile") << "Processing read_request: \n"
                                   << request->DebugString() << "\n";
  
  // Access the data in the index.  Return the data on success, or
  // notify the remote client of our failure if we did not succeed.
  try {
    IndexRecord record = m_index.read_request(request->col(), 
                                              request->row(), 
                                              request->depth(), 
                                              request->transaction_id());
    *(response->mutable_index_record()) = record;
  } catch (TileNotFoundErr &e) {
    controller->SetFailed(e.what());
  }
  
  done->Run();
}

