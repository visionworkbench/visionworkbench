// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/common.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/IndexService.h>

#include <google/protobuf/descriptor.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

void null_closure() {}

// -----------------------------------------------------------------------------
//                                IndexServer
// -----------------------------------------------------------------------------

class IndexServer : public google::protobuf::RpcController {
  std::string m_exchange;
  std::string m_queue;
  AmqpConnection m_conn;
  bool m_failed;
  std::string m_failed_reason;
  
  boost::shared_ptr<google::protobuf::Service> m_service;

public:
  IndexServer(std::string exchange, std::string queue) :
    m_exchange(exchange), m_queue(queue) {

    m_conn.exchange_declare(exchange, "direct", true, false);
    m_conn.queue_declare(queue, true, true, false);

    this->Reset();
  }

  virtual ~IndexServer() {}

  void export_with_routing_key(boost::shared_ptr<google::protobuf::Service> service,
                               std::string routing_key) {
    m_conn.queue_bind(m_queue, m_exchange, routing_key);
    std::cout << "\t    Bound \'" << m_queue << "\' to \'" << m_exchange 
              << "\' with the \'" << routing_key << "\' routing key.\n";
    m_service = service;
  }

  void run() {

    while(1) {

      try {

        // Step 1 : Wait for an incoming message.
        std::string routing_key;
        std::cout << "\n\nWaiting for message...\n";
        boost::shared_array<uint8> request_bytes = m_conn.basic_consume(m_queue, 
                                                                        routing_key, 
                                                                        false);
        WireMessage wire_request(request_bytes);
        RpcRequestWrapper request_wrapper = wire_request.parse_as_message<RpcRequestWrapper>();
        std::cout << "[RPC] " << request_wrapper.method() 
                  << " from " << request_wrapper.requestor() << "\n";
        
        // Step 2 : Delegate the message to the proper method on the service.
        const google::protobuf::MethodDescriptor* method = 
          m_service->GetDescriptor()->FindMethodByName(request_wrapper.method());
        
        boost::shared_ptr<google::protobuf::Message> 
          request(m_service->GetRequestPrototype(method).New());
        
        boost::shared_ptr<google::protobuf::Message> 
          response(m_service->GetResponsePrototype(method).New());
        
        request->ParseFromString(request_wrapper.payload());
        std::cout << "Request:\n" << request->DebugString() << "\n";
        m_service->CallMethod(method, this, request.get(), response.get(), 
                              google::protobuf::NewCallback(&null_closure));
        std::cout << "Response:\n" << response->DebugString() << "\n\n";
        
        // Step 3 : Return the result.
        RpcResponseWrapper response_wrapper;
        response_wrapper.set_payload(response->SerializeAsString());
        response_wrapper.set_error(false);
        if (this->Failed()) {
          response_wrapper.set_error(true);
          response_wrapper.set_error_string(this->ErrorText());
          this->Reset();
        }
        
        //        std::cout << "Response wrapper " << response_wrapper.DebugString() << "\n";
        
        WireMessage wire_response(&response_wrapper);
        m_conn.basic_publish(wire_response.serialized_bytes(), 
                             wire_response.size(),
                             m_exchange, request_wrapper.requestor() );
      

      } catch (vw::Exception &e) {
        std::cout << "WARNING: A server error occurred: " << e.what() << "\n";
      } 

    }    
  }

  // ------------------ RpcController Methods ---------------------

  // Resets the RpcController to its initial state so that it may be reused in
  // a new call.  Must not be called while an RPC is in progress.
  virtual void Reset() {
    m_failed = false;
    m_failed_reason = "";
  }

  // After a call has finished, returns true if the call failed.  The possible
  // reasons for failure depend on the RPC implementation.  Failed() must not
  // be called before a call has finished.  If Failed() returns true, the
  // contents of the response message are undefined.
  virtual bool Failed() const {
    return m_failed;
  }
  
  // If Failed() is true, returns a human-readable description of the error.
  virtual std::string ErrorText() const {
    return m_failed_reason;
  }

  // Advises the RPC system that the caller desires that the RPC call be
  // canceled.  The RPC system may cancel it immediately, may wait awhile and
  // then cancel it, or may not even cancel the call at all.  If the call is
  // canceled, the "done" callback will still be called and the RpcController
  // will indicate that the call failed at that time.
  virtual void StartCancel() {}  

  // Server-side methods ---------------------------------------------
  
  // Causes Failed() to return true on the client side.  "reason" will be
  // incorporated into the message returned by ErrorText().  If you find
  // you need to return machine-readable information about failures, you
  // should incorporate it into your response protocol buffer and should
  // NOT call SetFailed().
  virtual void SetFailed(const std::string& reason) {
    m_failed = true;
    m_failed_reason = reason;
  }

  // If true, indicates that the client canceled the RPC, so the server may
  // as well give up on replying to it.  The server should still call the
  // final "done" callback.
  virtual bool IsCanceled() const { return false; }
  
  // Asks that the given callback be called when the RPC is canceled.  The
  // callback will always be called exactly once.  If the RPC completes without
  // being canceled, the callback will be called after completion.  If the RPC
  // has already been canceled when NotifyOnCancel() is called, the callback
  // will be called immediately.
  //
  // NotifyOnCancel() must be called no more than once per request.
  virtual void NotifyOnCancel(google::protobuf::Closure* callback) {}

};


// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string queue_name, root_directory;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value(""), "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("root-directory", po::value<std::string>(&root_directory));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("root-directory", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-q <queue name>] root_directory" <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("root-directory") != 1 ) {
    std::cerr << "Error: must specify a root directory that contains plate files!" 
              << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }
  
  std::cout << "Initializing the index server.\n";
  IndexServer server(INDEX_EXCHANGE, INDEX_QUEUE);

  std::cout << "\t--> Exporting service.\n";
  boost::shared_ptr<google::protobuf::Service> service( new IndexServiceImpl(root_directory) );
  server.export_with_routing_key(service, "index");

  std::cout << "\t--> Listening for messages.\n";
  server.run();

  return 0;
}

