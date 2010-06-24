#ifndef __VW_PLATE_MOD_PLATE_HANDLERS_H__
#define __VW_PLATE_MOD_PLATE_HANDLERS_H__

struct request_rec;

namespace vw {
namespace platefile {
  class QueryMap;
}}

int handle_image(request_rec *r, const std::string& url, const vw::platefile::QueryMap& query);
int  handle_wtml(request_rec *r, const std::string& url, const vw::platefile::QueryMap& query);

#endif
