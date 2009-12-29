// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Image/PixelTypeInfo.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/LocalIndex.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/PagedIndex.h>
using namespace vw::platefile;
using namespace vw;

// -------------------------------------------------------------------
//                            PLATE INDEX
// -------------------------------------------------------------------

boost::shared_ptr<Index> Index::construct_create(std::string url, IndexHeader new_index_info) {

  // Remote index URLs start with "pf://..."
  if (url.find("pf://") == 0) {
    return boost::shared_ptr<Index>(new RemoteIndex(url, new_index_info));

  // All other URLs are passed along to the local index structure
  } else {
    //return boost::shared_ptr<Index>(new LocalTreeIndex(url, new_index_info));
    return boost::shared_ptr<Index>(new PagedIndex(url, new_index_info));
  }

}

boost::shared_ptr<Index> Index::construct_open(std::string url) {

  // Remote index URLs start with "pf://..."
  if (url.find("pf://") == 0) {
    return boost::shared_ptr<Index>(new RemoteIndex(url));

  // All other URLs are passed along to the local index structure
  } else {
    //return boost::shared_ptr<Index>(new LocalTreeIndex(url));
    return boost::shared_ptr<Index>(new PagedIndex(url));
  }

}

