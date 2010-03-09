// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/Exception.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/common.h>
#include <vw/Plate/ProtoBuffers.pb.h>
using namespace vw;
using namespace vw::platefile;

#include <unistd.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
 
// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// Parse a URL with the format: pf://<exchange>/<platefile name>.plate
void parse_url(std::string const& url, std::string &hostname, int &port, 
               std::string &exchange, std::string &platefile_name) {

    std::string bare_exchange;
    if (url.find("pf://") != 0) {
      vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- this does not appear to be a well-formed URL: " << url);
    } else {
      std::string substr = url.substr(5);

      std::vector<std::string> split_vec;
      boost::split( split_vec, substr, boost::is_any_of("/") );

      // No hostname was specified: pf://<exchange>/<platefilename>.plate
      if (split_vec.size() == 2) {
        hostname = "localhost"; // default to localhost
        port = 5672;            // default rabbitmq port

        bare_exchange = split_vec[0];
        platefile_name = split_vec[1];

      // No hostname was specified: pf://<ip address>:<port>/<exchange>/<platefilename>.plate
      } else if (split_vec.size() == 3) {

        bare_exchange = split_vec[1];
        platefile_name = split_vec[2];

        std::vector<std::string> host_port_vec;
        boost::split( host_port_vec, split_vec[0], boost::is_any_of(":") );

        if (host_port_vec.size() == 1) {

          hostname = host_port_vec[0];
          port = 5672;            // default rabbitmq port

        } else if (host_port_vec.size() == 2) {

          hostname = host_port_vec[0];
          port = boost::lexical_cast<int>(host_port_vec[1]);

        } else {
          vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- " 
                   << "could not parse hostname and port from URL string: " << url);
        }

      } else {
        vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- "
                 << "could not parse URL string: " << url);
      }
    }

    // From here, we'll always be operating within the platefile exchange
    // namespace... so prepend it.
    exchange = std::string(PLATE_EXCHANGE_NAMESPACE) + "." + bare_exchange;
}

// ----------------------------------------------------------------------
//                         REMOTE INDEX PAGE
// ----------------------------------------------------------------------

vw::platefile::RemoteIndexPage::RemoteIndexPage(int platefile_id, 
                                                boost::shared_ptr<AmqpRpcClient> rpc_controller,
                                                boost::shared_ptr<IndexService> index_service,
                                                int level, int base_col, int base_row, 
                                                int page_width, int page_height) :
  IndexPage(level, base_col, base_row, page_width, page_height),
  m_platefile_id(platefile_id), m_rpc_controller(rpc_controller),
  m_index_service(index_service) {

  // Use the PageRequest RPC to fetch the remote page from the index
  // server.
  IndexPageRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_col(base_col);
  request.set_row(base_row);
  request.set_level(level);
  IndexPageReply response;
  try {
    m_index_service->PageRequest(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));
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

vw::platefile::RemoteIndexPage::~RemoteIndexPage() {
  this->sync();
}

// Hijack this method momentartarily to mark the page as "dirty" by
// setting m_needs_saving to true.
void vw::platefile::RemoteIndexPage::set(TileHeader const& header, IndexRecord const& record) {
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

void vw::platefile::RemoteIndexPage::flush_write_queue() {
  if (m_write_queue.size() > 0) {

    // For debugging:
    //    std::cout << "Call to the new flush_write_queue() with " 
    //              << m_write_queue.size() << " entries.\n";

    IndexMultiWriteUpdate request;
    while (m_write_queue.size() > 0) {
      *(request.mutable_write_updates()->Add()) = m_write_queue.front();
      m_write_queue.pop();
    }
    RpcNullMessage response;
    m_index_service->MultiWriteUpdate(m_rpc_controller.get(), &request, &response, 
                                      google::protobuf::NewCallback(&null_closure));

  }
}

void vw::platefile::RemoteIndexPage::sync() {
  this->flush_write_queue();
}

// ----------------------------------------------------------------------
//                    REMOTE INDEX PAGE GENERATOR
// ----------------------------------------------------------------------

vw::platefile::RemotePageGenerator::RemotePageGenerator( int platefile_id, 
                                                         boost::shared_ptr<AmqpRpcClient> rpc_controller,
                                                         boost::shared_ptr<IndexService> index_service,
                                                         int level, int base_col, int base_row, 
                                                         int page_width, int page_height) : 
  m_platefile_id(platefile_id), m_rpc_controller(rpc_controller), 
  m_index_service(index_service), m_level(level), 
  m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height) {}  

boost::shared_ptr<vw::platefile::IndexPage> 
vw::platefile::RemotePageGenerator::generate() const {
  return boost::shared_ptr<IndexPage>(new RemoteIndexPage(m_platefile_id, m_rpc_controller,
                                                          m_index_service, 
                                                          m_level, m_base_col, m_base_row,
                                                          m_page_width, m_page_height) );
}

boost::shared_ptr<IndexPageGenerator> RemotePageGeneratorFactory::create(int level, 
                                                                         int base_col, 
                                                                         int base_row, 
                                                                         int page_width, 
                                                                         int page_height) {
  VW_ASSERT(m_platefile_id != -1 && m_rpc_controller && m_index_service,
            LogicErr() << "Error: RemotePageGeneratorFactory has not yet been initialized.");

  // Create the proper type of page generator.
  boost::shared_ptr<RemotePageGenerator> page_gen;
  page_gen.reset( new RemotePageGenerator(m_platefile_id, m_rpc_controller,
                                          m_index_service,
                                          level, base_col, base_row,
                                          page_width, page_height) );
  
  // Wrap it in tho IndexPageGenerator class and return it.
  return boost::shared_ptr<IndexPageGenerator>( new IndexPageGenerator(page_gen) );
}

// ----------------------------------------------------------------------
//                             REMOTE INDEX
// ----------------------------------------------------------------------


/// Constructor (for opening an existing Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url) :
  PagedIndex(boost::shared_ptr<PageGeneratorFactory>( new RemotePageGeneratorFactory() ))  {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string hostname;
  int port;
  std::string exchange;
  std::string platefile_name;
  parse_url(url, hostname, port, exchange, platefile_name);

  std::string queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection(hostname, port));
  m_rpc_controller.reset( new AmqpRpcClient(conn, exchange, queue_name, "index") );
  m_index_service.reset ( new IndexService::Stub(m_rpc_controller.get() ) );
  m_rpc_controller->bind_service(m_index_service, queue_name);

  // Send an IndexOpenRequest to the AMQP index server.
  IndexOpenRequest request;
  request.set_plate_name(platefile_name);

  IndexOpenReply response;
  m_index_service->OpenRequest(m_rpc_controller.get(), &request, &response,
                               google::protobuf::NewCallback(&null_closure));

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();

  // Properly initialize the PageGenFactory and set it.
  boost::shared_ptr<PageGeneratorFactory> factory( new RemotePageGeneratorFactory(m_platefile_id,
                                                                                  m_rpc_controller, 
                                                                                  m_index_service));
  this->set_page_generator_factory(factory);

  // Every time you run num_levels, it synchronizes the number of
  // local (cached) levels with the number of levels on the index
  // server, so we run it here to synchronize it for the first time.
  this->num_levels();

  vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile_name
                               << "\"   ID: " << m_platefile_id << "\n";

  vw_out() << "Opened remote platefile \"" << platefile_name
            << "\"   ID: " << m_platefile_id << "  ( " << this->num_levels() << " levels )\n";
}

/// Constructor (for creating a new Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url, IndexHeader index_header_info) :
  PagedIndex(boost::shared_ptr<PageGeneratorFactory>( new RemotePageGeneratorFactory() )),
  m_index_header(index_header_info) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string hostname;
  int port;
  std::string exchange;
  std::string platefile_name;
  parse_url(url, hostname, port, exchange, platefile_name);

  std::string queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection(hostname, port));
  m_rpc_controller.reset( new AmqpRpcClient(conn, exchange, queue_name, "index") );
  m_index_service.reset ( new IndexService::Stub(m_rpc_controller.get() ) );
  m_rpc_controller->bind_service(m_index_service, queue_name);

  // Send an IndexCreateRequest to the AMQP index server.
  IndexCreateRequest request;
  request.set_plate_name(platefile_name);
  index_header_info.set_platefile_id(0);  // this takes care of this 'required' 
                                          // protobuf property, which is not set yet
  *(request.mutable_index_header()) = index_header_info;

  IndexOpenReply response;
  m_index_service->CreateRequest(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();

  // Properly initialize the PageGenFactory and set it.
  boost::shared_ptr<PageGeneratorFactory> factory( new RemotePageGeneratorFactory(m_platefile_id,
                                                                                  m_rpc_controller, 
                                                                                  m_index_service));
  this->set_page_generator_factory(factory);

  // Every time you run num_levels, it synchronizes the number of
  // local (cached) levels with the number of levels on the index
  // server, so we run it here to synchronize it for the first time.
  this->num_levels();

  vw_out(InfoMessage, "plate") << "Created remote platefile \"" << platefile_name
                               << "\"   ID: " << m_platefile_id << "\n";
  vw_out() << "Created remote platefile \"" << platefile_name
           << "\"   ID: " << m_platefile_id << "  ( " << this->num_levels() << " levels )\n";

}


/// Destructor
vw::platefile::RemoteIndex::~RemoteIndex() {}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::RemoteIndex::write_request(uint64 &size) {
  IndexWriteRequest request;
  request.set_platefile_id(m_platefile_id);

  IndexWriteReply response;
  m_index_service->WriteRequest(m_rpc_controller.get(), &request, &response, 
                                google::protobuf::NewCallback(&null_closure));
  size = response.size();
  return response.blob_id();
}

/// Log a message to the platefile log.
void vw::platefile::RemoteIndex::log(std::string message) {
  IndexLogRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_message(message);

  RpcNullMessage response;
  m_index_service->LogRequest(m_rpc_controller.get(), &request, &response, 
                       google::protobuf::NewCallback(&null_closure));  
}

/// Writing, pt. 3: Signal the completion 
void vw::platefile::RemoteIndex::write_complete(int blob_id, uint64 blob_offset) {

  // First we make sure that we flush the write queue by synchronizing
  // all of the pages back to the index_server!  Otherwise the write
  // won't actually be complete!!
  this->sync();

  // Then we issue an actual write_complete.
  IndexWriteComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_blob_id(blob_id);
  request.set_blob_offset(blob_offset);

  RpcNullMessage response;
  m_index_service->WriteComplete(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));
}

// This is the old valid_tiles RPC implementation.  It has been
// replaced by the valid_tiles implementation in the PagedIndex class.
// 
// std::list<TileHeader> vw::platefile::RemoteIndex::valid_tiles(int level, BBox2i const& region,
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
//   m_index_service->ValidTiles(m_rpc_controller.get(), &request, &response, 
//                               google::protobuf::NewCallback(&null_closure));

//   std::list<TileHeader> results;
//   for (int i = 0; i < response.tile_headers_size(); ++i) {
//     results.push_back(response.tile_headers().Get(i));    
//   }
//   return results;
// }

vw::int32 vw::platefile::RemoteIndex::num_levels() const { 
  IndexNumLevelsRequest request;
  request.set_platefile_id(m_platefile_id);
  IndexNumLevelsReply response;
  m_index_service->NumLevelsRequest(m_rpc_controller.get(), &request, &response, 
                                    google::protobuf::NewCallback(&null_closure));

  // Make sure that the local (cached) number of levels matches the
  // number of levels on the server.
  for (int level = m_levels.size(); level < response.num_levels(); ++level) { 
    boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_page_gen_factory, level, 
                                                            m_page_width, m_page_height, 
                                                            m_default_cache_size) );
    m_levels.push_back(new_level);
  }

  // Return the number of levels.
  return response.num_levels();
}

vw::int32 vw::platefile::RemoteIndex::version() const { 
  return m_index_header.version();
}

std::string vw::platefile::RemoteIndex::platefile_name() const { 
  return m_full_plate_filename;
}

IndexHeader RemoteIndex::index_header() const { 
  return m_index_header;
}

vw::int32 vw::platefile::RemoteIndex::tile_size() const { 
  return m_index_header.tile_size();
}

std::string vw::platefile::RemoteIndex::tile_filetype() const { 
  return m_index_header.tile_filetype();
}

vw::PixelFormatEnum vw::platefile::RemoteIndex::pixel_format() const {
  return vw::PixelFormatEnum(m_index_header.pixel_format());
}

vw::ChannelTypeEnum vw::platefile::RemoteIndex::channel_type() const {
  return vw::ChannelTypeEnum(m_index_header.channel_type());
}

// --------------------- TRANSACTIONS ------------------------

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
vw::int32 vw::platefile::RemoteIndex::transaction_request(std::string transaction_description,
                                                          int transaction_id_override) {

  IndexTransactionRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_description(transaction_description);
  request.set_transaction_id_override(transaction_id_override);
  
  IndexTransactionReply response;
  m_index_service->TransactionRequest(m_rpc_controller.get(), &request, &response, 
                                      google::protobuf::NewCallback(&null_closure));
  return response.transaction_id();
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuding a transaction_complete method.
void vw::platefile::RemoteIndex::transaction_complete(int32 transaction_id, bool update_read_cursor) {
  IndexTransactionComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  request.set_update_read_cursor(update_read_cursor);
  
  RpcNullMessage response;
  m_index_service->TransactionComplete(m_rpc_controller.get(), &request, &response, 
                                       google::protobuf::NewCallback(&null_closure));
}

// If a transaction fails, we may need to clean up the mosaic.  
void vw::platefile::RemoteIndex::transaction_failed(int32 transaction_id) {
  IndexTransactionFailed request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  
  RpcNullMessage response;
  m_index_service->TransactionFailed(m_rpc_controller.get(), &request, &response, 
                                     google::protobuf::NewCallback(&null_closure));
}

vw::int32 vw::platefile::RemoteIndex::transaction_cursor() {
  IndexTransactionCursorRequest request;
  request.set_platefile_id(m_platefile_id);
  
  IndexTransactionCursorReply response;
  m_index_service->TransactionCursor(m_rpc_controller.get(), &request, &response, 
                                     google::protobuf::NewCallback(&null_closure));
  return response.transaction_id();
}


