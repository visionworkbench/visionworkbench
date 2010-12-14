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
#include <boost/foreach.hpp>
namespace fs = boost::filesystem;
namespace pb = ::google::protobuf;

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

IndexServiceImpl::IndexServiceRecord* IndexServiceImpl::add_index(std::string root_directory,
                                                                  std::string plate_filename,
                                                                  boost::shared_ptr<Index> index) {

    // Build up an IndexServiceRecord
    IndexServiceRecord rec;
    rec.short_plate_filename = plate_filename;
    rec.full_plate_filename = root_directory + "/" + plate_filename;
    rec.index = index;

    // Store the record in a std::map by platefile_id
    m_indices[index->index_header().platefile_id()]  = rec;
    return &m_indices[index->index_header().platefile_id()];
}


IndexServiceImpl::IndexServiceRecord* IndexServiceImpl::find_id(int32 platefile_id) {
  index_list_type::iterator i = m_indices.find(platefile_id);
  if (i == m_indices.end())
    return 0;
  return &i->second;
}

IndexServiceImpl::IndexServiceRecord IndexServiceImpl::find_id_throw(int32 platefile_id) {
  IndexServiceRecord *rec = find_id(platefile_id);
  if (rec)
    return *rec;
  vw_throw(InvalidPlatefileErr() << "No platefile matching this platefile id found: " << platefile_id);
}

IndexServiceImpl::IndexServiceRecord* IndexServiceImpl::find_name(const std::string& name) {
  BOOST_FOREACH(index_list_type::value_type& index, m_indices)
    if ( name == index.second.short_plate_filename )
      return &(index.second);
  return 0;
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
  BOOST_FOREACH(index_list_type::value_type& i, m_indices) {
    vw_out() << "\t--> Syncing index for " << i.second.short_plate_filename << " to disk.\n";
    i.second.index->sync();
  }
}

#define METHOD_IMPL(Name, Input, Output) \
  void IndexServiceImpl::Name(pb::RpcController*, const Input* request, Output* response, pb::Closure* done)
#define METHOD_IMPL_NOREPLY(Name, Input) \
  void IndexServiceImpl::Name(pb::RpcController*, const Input* request, RpcNullMsg*, pb::Closure* done)

METHOD_IMPL(OpenRequest, IndexOpenRequest, IndexOpenReply) {
  detail::RequireCall call(done);
  IndexServiceRecord *r = find_name(request->plate_name());

  if (!r)
    vw_throw(InvalidPlatefileErr() << "No platefile matching this platename found: " << request->plate_name());

  response->set_short_plate_filename( r->short_plate_filename );
  response->set_full_plate_filename( r->full_plate_filename );
  *(response->mutable_index_header()) = r->index->index_header();
}

METHOD_IMPL(CreateRequest, IndexCreateRequest, IndexOpenReply) {
  detail::RequireCall call(done);
  IndexServiceRecord *r = find_name(request->plate_name());

  const Url url(m_root_directory + "/" + request->plate_name());
  try {
    if (!r && !fs::exists(url.path())) {
      // CREATE: If there was no index already open by that name, and
      // no platefile on disk, then we try to create it.
      fs::create_directory(url.path());
      boost::shared_ptr<Index> idx = Index::construct_create(url, request->index_header());
      r = this->add_index(m_root_directory, request->plate_name(), idx);
    } else {
      if (!r) {
      // REOPEN: If a platefile actually does appear to exist on disk,
      // but we had not previously opened it, we open it here.
      boost::shared_ptr<Index> idx = Index::construct_open(url);
      r = this->add_index(m_root_directory, request->plate_name(), idx);
    }
    // OPEN: Now check that the opened/reopened index matches the requested format.
    const IndexHeader& orig = r->index->index_header();
    const IndexHeader& req  = request->index_header();

    // Most of these settings are just defaults or can be gleaned from the
    // TileHeader or the data itself (like, the tile size). However, the type
    // cannot be changed.
#define CHECK(field) VW_ASSERT(orig.field() == req.field(), PlatefileCreationErr() << "Refusing to override " #field " on platefile " << request->plate_name())
    CHECK(type);
#undef CHECK
    }
  } catch (const IOErr& e) {
    vw_throw(PlatefileCreationErr() << e.what());
  }

  response->set_short_plate_filename( r->short_plate_filename );
  response->set_full_plate_filename( r->full_plate_filename );
  *(response->mutable_index_header()) = r->index->index_header();
}

METHOD_IMPL(InfoRequest, IndexInfoRequest, IndexInfoReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());

  response->set_short_plate_filename(rec.short_plate_filename);
  response->set_full_plate_filename(rec.full_plate_filename);
  *(response->mutable_index_header()) = rec.index->index_header();
}

METHOD_IMPL(ListRequest, IndexListRequest, IndexListReply) {
  detail::RequireCall call(done);
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
}

METHOD_IMPL(PageRequest, IndexPageRequest, IndexPageReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());

  boost::shared_ptr<IndexPage> page =
    rec.index->page_request(request->col(), request->row(), request->level());

  std::ostringstream ostr;
  page->serialize(ostr);
  response->set_page_bytes(ostr.str().c_str(), ostr.str().size());
}

METHOD_IMPL(ReadRequest, IndexReadRequest, IndexReadReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());

  *(response->mutable_index_record()) =
    rec.index->read_request(request->col(), request->row(), request->level(),
                            request->transaction_id(), request->exact_transaction_match());
}

METHOD_IMPL(WriteRequest, IndexWriteRequest, IndexWriteReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());

  uint64 size;
  int blob_id = rec.index->write_request(size);
  response->set_blob_id(blob_id);
  response->set_size(size);
}

METHOD_IMPL_NOREPLY(WriteUpdate, IndexWriteUpdate) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  rec.index->write_update(request->header(), request->record());
}

METHOD_IMPL_NOREPLY(MultiWriteUpdate, IndexMultiWriteUpdate) {
  detail::RequireCall call(done);
  // MultiWrite updates are packetized.  We iterate over them here.
  for (int i = 0; i < request->write_updates().size(); ++i) {
    IndexWriteUpdate update = request->write_updates().Get(i);
    IndexServiceRecord rec = find_id_throw(update.platefile_id());
    rec.index->write_update(update.header(), update.record());
  }
}

METHOD_IMPL_NOREPLY(WriteComplete, IndexWriteComplete) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  rec.index->write_complete(request->blob_id(), request->blob_offset());
}

METHOD_IMPL(TransactionRequest, IndexTransactionRequest, IndexTransactionReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  Transaction transaction_id = 
    rec.index->transaction_request(request->description(), request->transaction_id_override());
  response->set_transaction_id(transaction_id);
}

METHOD_IMPL_NOREPLY(TransactionComplete, IndexTransactionComplete) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  rec.index->transaction_complete(request->transaction_id(), request->update_read_cursor());
}

METHOD_IMPL_NOREPLY(TransactionFailed, IndexTransactionFailed) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  rec.index->transaction_failed(request->transaction_id());
}

METHOD_IMPL(TransactionCursor, IndexTransactionCursorRequest, IndexTransactionCursorReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  Transaction transaction_id = rec.index->transaction_cursor();
  response->set_transaction_id(transaction_id);
}

METHOD_IMPL(NumLevelsRequest, IndexNumLevelsRequest, IndexNumLevelsReply) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  response->set_num_levels(rec.index->num_levels());
}

METHOD_IMPL_NOREPLY(LogRequest, IndexLogRequest) {
  detail::RequireCall call(done);
  IndexServiceRecord rec = find_id_throw(request->platefile_id());
  rec.index->log(request->message());
}

METHOD_IMPL(TestRequest, IndexTestRequest, IndexTestReply) {
  detail::RequireCall call(done);
  response->set_value(request->value());
}

