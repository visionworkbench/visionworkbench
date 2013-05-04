// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Plate/Datastore.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/detail/Blobstore.h>
#include <vw/Plate/detail/Dirstore.h>

namespace vw { namespace platefile {

Datastore* Datastore::open(const Url& url) {
  if (url.scheme() == "dir")
    return new detail::Dirstore(url);
  else
    return new detail::Blobstore(url);
  //vw_throw(NoImplErr() << "Unsupported datastore scheme: " << url.scheme());
}

Datastore* Datastore::open(const Url& url, const IndexHeader& d) {
  if (url.scheme() == "dir")
    return new detail::Dirstore(url, d);
  else
    return new detail::Blobstore(url, d);
  //vw_throw(NoImplErr() << "Unsupported datastore scheme: " << url.scheme());
}

Datastore::TileSearch& Datastore::get(TileSearch& buf, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit) {
  this->head(buf, level, row, col, range, limit);
  return this->populate(buf);
}

Datastore::TileSearch& Datastore::get(TileSearch& buf, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit) {
  this->head(buf, level, region, range, limit);
  return this->populate(buf);
}

}}
