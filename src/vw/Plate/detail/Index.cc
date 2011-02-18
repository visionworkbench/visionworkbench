// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

