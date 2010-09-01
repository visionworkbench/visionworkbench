// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


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
