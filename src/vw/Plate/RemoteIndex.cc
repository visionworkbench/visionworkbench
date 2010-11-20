// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/Rpc.h>
#include <vw/Plate/IndexService.pb.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/Exception.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <unistd.h>

std::string split_url(Url& url) {
  Url::split_t sp = url.path_split();
  VW_ASSERT(sp.size() > 0, ArgumentErr() << "Expected a platefile url (bad path)");

  const std::string platefile_name = sp.back();
  VW_ASSERT(boost::ends_with(platefile_name, ".plate"), ArgumentErr() << "Expected a platefile url (doesn't end in .plate)");
  sp.pop_back();
  url.path_join(sp);
  return platefile_name;
}

RemoteIndexPage::RemoteIndexPage(int platefile_id,
                                 boost::shared_ptr<IndexClient> client,
                                 int level, int base_col, int base_row,
                                 int page_width, int page_height)
  : IndexPage(level, base_col, base_row, page_width, page_height),
    m_platefile_id(platefile_id), m_client(client)
{
  // Use the PageRequest RPC to fetch the remote page from the index
  // server.
  IndexPageRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_col(base_col);
  request.set_row(base_row);
  request.set_level(level);
  IndexPageReply response;
  try {
    m_client->PageRequest(m_client.get(), &request, &response, null_callback());
  } catch (TileNotFoundErr &e) {
    // If the remote server doesn't have this index page, that's ok.
    // Doing nothing here will create an empty index page in this
    // local cache.  If we end up writing data to the local page, then
    // the page will be created once we start writing some of our data
    // to the index server as needed.
  }

  // Deserialize the data from the message to populate this index page.
  std::istringstream istr(response.page_bytes(), std::ios::binary);
  this->deserialize(istr);
}

RemoteIndexPage::~RemoteIndexPage() {
  this->sync();
}

// Hijack this method momentartarily to mark the page as "dirty" by
// setting m_needs_saving to true.
void RemoteIndexPage::set(TileHeader const& header, IndexRecord const& record) {
  // First call up to the parent class and let the original code run.
  IndexPage::set(header, record);

  // Save this write request to the queue, and flush the queue if it
  // has gotten too full.  (We send an update to the index_server
  // every 50 writes!)
  IndexWriteUpdate request;
  request.set_platefile_id(m_platefile_id);
  *(request.mutable_header()) = header;
  *(request.mutable_record()) = record;
  m_write_queue.push(request);
  static const unsigned write_queue_size = 50; // arbitrary, but probably good.
  if (m_write_queue.size() >= write_queue_size)
    this->flush_write_queue();
}

void RemoteIndexPage::flush_write_queue() {
  if (m_write_queue.size() > 0) {

    // For debugging:
    //    std::cout << "Call to the new flush_write_queue() with "
    //              << m_write_queue.size() << " entries.\n";

    IndexMultiWriteUpdate request;
    while (m_write_queue.size() > 0) {
      *(request.mutable_write_updates()->Add()) = m_write_queue.front();
      m_write_queue.pop();
    }
    RpcNullMsg response;
    m_client->MultiWriteUpdate(m_client.get(), &request, &response, null_callback());

  }
}

void RemoteIndexPage::sync() {
  this->flush_write_queue();
}

// ----------------------------------------------------------------------
//                    REMOTE INDEX PAGE GENERATOR
// ----------------------------------------------------------------------

RemotePageGenerator::RemotePageGenerator( int platefile_id,
                                          boost::shared_ptr<IndexClient> client,
                                          int level, int base_col, int base_row,
                                          int page_width, int page_height)
  : m_platefile_id(platefile_id), m_client(client), m_level(level),
    m_base_col(base_col), m_base_row(base_row),
    m_page_width(page_width), m_page_height(page_height) {}

boost::shared_ptr<IndexPage>
RemotePageGenerator::generate() const {
  return boost::shared_ptr<IndexPage>(
      new RemoteIndexPage(m_platefile_id, m_client, m_level,
                          m_base_col, m_base_row, m_page_width, m_page_height) );
}

boost::shared_ptr<PageGeneratorBase>
RemotePageGeneratorFactory::create(int level, int base_col, int base_row, int page_width, int page_height) {
  //VW_ASSERT(m_platefile_id != -1 && m_rpc_controller && m_index_service,
  //          LogicErr() << "Error: RemotePageGeneratorFactory has not yet been initialized.");

  // Create the proper type of page generator.
  boost::shared_ptr<PageGeneratorBase> page_gen(
    new RemotePageGenerator(m_platefile_id, m_client,
                            level, base_col, base_row,
                            page_width, page_height) );

  return page_gen;
}

std::string RemotePageGeneratorFactory::who() const {
  VW_ASSERT(m_platefile_id != -1,
            LogicErr() << "Error: RemotePageGeneratorFactory has not yet been initialized.");
  return vw::stringify(m_platefile_id);
}

// ----------------------------------------------------------------------
//                             REMOTE INDEX
// ----------------------------------------------------------------------


// Constructor (for opening an existing Index)
// expecting a url in the form of scheme://hostname:port/path/to/server/platefile.plate
RemoteIndex::RemoteIndex(const Url& url_)
  : m_url(url_), m_short_plate_filename(split_url(m_url)),
    m_client(new IndexClient(m_url))
{
  // TODO: some way to set client_name from here?
  // m_client->conn(url, string("remote_index_") + platefile_name);

  IndexOpenRequest request;
  request.set_plate_name(m_short_plate_filename);

  IndexOpenReply response;
  m_client->OpenRequest(m_client.get(), &request, &response, null_callback());

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();

  // Properly initialize the PageGenFactory and set it.
  boost::shared_ptr<PageGeneratorFactory> factory(
      new RemotePageGeneratorFactory(m_platefile_id, m_client));

  this->set_page_generator_factory(factory);
  this->set_default_cache_size(m_url.query().get("cache_size", 100u));

  // Every time you run num_levels, it synchronizes the number of
  // local (cached) levels with the number of levels on the index
  // server, so we run it here to synchronize it for the first time.
  this->num_levels();

  vw_out(InfoMessage, "plate")
    << "Opened remote platefile name[" << m_short_plate_filename
    << "] id[" << m_platefile_id << "]" << std::endl;
}

/// Constructor (for creating a new Index)
RemoteIndex::RemoteIndex(const Url& url_, IndexHeader index_header_info)
  : m_url(url_), m_short_plate_filename(split_url(m_url)),
    m_index_header(index_header_info),
    m_client(new IndexClient(m_url))
{
  // TODO: some way to set client_name from here?
  // m_client->conn(url, string("remote_index_") + platefile_name);

  // Send an IndexCreateRequest to the AMQP index server.
  IndexCreateRequest request;
  request.set_plate_name(m_short_plate_filename);
  index_header_info.set_platefile_id(0);  // this takes care of this 'required'
                                          // protobuf property, which is not set yet
  *(request.mutable_index_header()) = index_header_info;

  IndexOpenReply response;
  m_client->CreateRequest(m_client.get(), &request, &response, null_callback());

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();

  // Properly initialize the PageGenFactory and set it.
  boost::shared_ptr<PageGeneratorFactory> factory(
      new RemotePageGeneratorFactory(m_platefile_id, m_client));

  this->set_page_generator_factory(factory);
  this->set_default_cache_size(m_url.query().get("cache_size", 100u));

  // Every time you run num_levels, it synchronizes the number of
  // local (cached) levels with the number of levels on the index
  // server, so we run it here to synchronize it for the first time.
  this->num_levels();

  vw_out(InfoMessage, "plate")
    << "Created remote platefile name[" << m_short_plate_filename
    << "] id[" << m_platefile_id << "] levels[" << this->num_levels() << "]" << std::endl;
}


/// Destructor
RemoteIndex::~RemoteIndex() {}

// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int RemoteIndex::write_request(uint64 &size) {
  IndexWriteRequest request;
  request.set_platefile_id(m_platefile_id);

  IndexWriteReply response;
  m_client->WriteRequest(m_client.get(), &request, &response, null_callback());
  size = response.size();
  return response.blob_id();
}

/// Log a message to the platefile log.
void RemoteIndex::log(std::string message) {
  IndexLogRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_message(message);

  RpcNullMsg response;
  m_client->LogRequest(m_client.get(), &request, &response, null_callback());
}

/// Writing, pt. 3: Signal the completion
void RemoteIndex::write_complete(int blob_id, uint64 blob_offset) {

  // First we make sure that we flush the write queue by synchronizing
  // all of the pages back to the index_server!  Otherwise the write
  // won't actually be complete!!
  this->sync();

  // Then we issue an actual write_complete.
  IndexWriteComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_blob_id(blob_id);
  request.set_blob_offset(blob_offset);

  RpcNullMsg response;
  m_client->WriteComplete(m_client.get(), &request, &response, null_callback());
}

// This is the old valid_tiles RPC implementation.  It has been
// replaced by the valid_tiles implementation in the PagedIndex class.
//
// std::list<TileHeader> RemoteIndex::valid_tiles(int level, BBox2i const& region,
//                                                               int begin_transaction_id,
//                                                               int end_transaction_id,
//                                                               int min_num_matches) const {
//   IndexValidTilesRequest request;
//   request.set_platefile_id(m_platefile_id);
//   request.set_level(level);
//   request.set_region_col(region.min().x());
//   request.set_region_row(region.min().y());
//   request.set_region_width(region.width());
//   request.set_region_height(region.height());
//   request.set_begin_transaction_id(begin_transaction_id);
//   request.set_end_transaction_id(end_transaction_id);
//   request.set_min_num_matches(min_num_matches);

//   IndexValidTilesReply response;
//   m_client->ValidTiles(m_client.get(), &request, &response, null_callback());

//   std::list<TileHeader> results;
//   for (int i = 0; i < response.tile_headers_size(); ++i) {
//     results.push_back(response.tile_headers().Get(i));
//   }
//   return results;
// }

vw::int32 RemoteIndex::num_levels() const {
  IndexNumLevelsRequest request;
  request.set_platefile_id(m_platefile_id);
  IndexNumLevelsReply response;
  m_client->NumLevelsRequest(m_client.get(), &request, &response, null_callback());

  // Make sure that the local (cached) number of levels matches the
  // number of levels on the server.
  for (int level = boost::numeric_cast<int>(m_levels.size()); level < response.num_levels(); ++level) {
    boost::shared_ptr<IndexLevel> new_level(
        new IndexLevel(m_page_gen_factory, level, m_page_width, m_page_height, m_default_cache_size) );
    m_levels.push_back(new_level);
  }

  // Return the number of levels.
  return response.num_levels();
}

vw::int32 RemoteIndex::version() const {
  return m_index_header.version();
}

std::string RemoteIndex::platefile_name() const {
  return m_full_plate_filename;
}

IndexHeader RemoteIndex::index_header() const {
  return m_index_header;
}

vw::int32 RemoteIndex::tile_size() const {
  return m_index_header.tile_size();
}

std::string RemoteIndex::tile_filetype() const {
  return m_index_header.tile_filetype();
}

vw::PixelFormatEnum RemoteIndex::pixel_format() const {
  return vw::PixelFormatEnum(m_index_header.pixel_format());
}

vw::ChannelTypeEnum RemoteIndex::channel_type() const {
  return vw::ChannelTypeEnum(m_index_header.channel_type());
}

// --------------------- TRANSACTIONS ------------------------

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
Transaction RemoteIndex::transaction_request(std::string transaction_description,
                                             TransactionOrNeg transaction_id_override) {
  int32 id;
  if (transaction_id_override.newest())
    id = -1;
  else
    id = static_cast<int32>(uint32(transaction_id_override.promote()));

  IndexTransactionRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_description(transaction_description);
  request.set_transaction_id_override(id);

  IndexTransactionReply response;
  m_client->TransactionRequest(m_client.get(), &request, &response, null_callback());
  return response.transaction_id();
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuding a transaction_complete method.
void RemoteIndex::transaction_complete(Transaction transaction_id, bool update_read_cursor) {
  IndexTransactionComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  request.set_update_read_cursor(update_read_cursor);

  RpcNullMsg response;
  m_client->TransactionComplete(m_client.get(), &request, &response, null_callback());
}

// If a transaction fails, we may need to clean up the mosaic.
void RemoteIndex::transaction_failed(Transaction transaction_id) {
  IndexTransactionFailed request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);

  RpcNullMsg response;
  m_client->TransactionFailed(m_client.get(), &request, &response, null_callback());
}

Transaction RemoteIndex::transaction_cursor() {
  IndexTransactionCursorRequest request;
  request.set_platefile_id(m_platefile_id);

  IndexTransactionCursorReply response;
  m_client->TransactionCursor(m_client.get(), &request, &response, null_callback());
  return response.transaction_id();
}


