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


#include <vw/Plate/detail/Index.h>
#include <vw/Plate/detail/RemoteIndex.h>
#include <vw/Plate/detail/LocalIndex.h>
#include <vw/Plate/HTTPUtils.h>

using namespace vw::platefile;
using namespace vw::platefile::detail;

boost::shared_ptr<Index> Index::construct_create(const Url& url, const IndexHeader& new_index_info) {
  if (url.scheme() == "file")
    return boost::shared_ptr<Index>(new LocalIndex(url.path(), new_index_info));
  else
    return boost::shared_ptr<Index>(new RemoteIndex(url, new_index_info));
}

boost::shared_ptr<Index> Index::construct_open(const Url& url) {
  if (url.scheme() == "file")
    return boost::shared_ptr<Index>(new LocalIndex(url.path()));
  else
    return boost::shared_ptr<Index>(new RemoteIndex(url));
}

