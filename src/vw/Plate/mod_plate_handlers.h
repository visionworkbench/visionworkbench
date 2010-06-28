#ifndef __VW_PLATE_MOD_PLATE_HANDLERS_H__
#define __VW_PLATE_MOD_PLATE_HANDLERS_H__

struct request_rec;

namespace vw {
namespace platefile {

class QueryMap;
class ApacheRequest;

int handle_image(const ApacheRequest& r);
int  handle_wtml(const ApacheRequest& r);

}} // namespace vw::platefile

#endif
