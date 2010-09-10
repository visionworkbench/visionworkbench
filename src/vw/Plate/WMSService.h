// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_INDEX_SERVICE_H__
#define __VW_PLATE_INDEX_SERVICE_H__

#include <vw/Plate/WMSMessages.pb.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/shared_ptr.hpp>

namespace vw {
namespace platefile {

  class PlateFile;

  class WMSServiceImpl : public WMSService {

    struct WMSServiceRecord {
      std::string plate_name;
      std::string full_path;
      boost::shared_ptr<PlateFile> platefile;
    };

    std::string m_root_directory;
    std::string m_cache_directory;

    typedef std::map<int32, WMSServiceRecord> platefile_list_type;
    platefile_list_type m_indices;

    // Private methods
    WMSServiceRecord add_platefile(std::string root_directory, std::string plate_filename,
                               boost::shared_ptr<PlateFile> platefile);

    /// Fetch an WMSServiceRecord for a given platefile_id, or throw an
    /// exception if no record is found.
    WMSServiceRecord get_platefile(int platefile_id);

    std::string create_image(const WMSTileRequest* req);

  public:

    WMSServiceImpl(std::string root_directory, std::string cache_directory);

    virtual void GetTile(::google::protobuf::RpcController* controller,
                         const ::vw::platefile::WMSTileRequest* request,
                         ::vw::platefile::WMSTileResponse* response,
                         ::google::protobuf::Closure* done);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_INDEX_SERVICE_H__
