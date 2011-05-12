// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/tools/Common.h>
#include <vw/FileIO/DiskImageResource.h>
#include <boost/scoped_ptr.hpp>

vw::ImageFormat vw::tools::taste_image(const std::string& filename) {
  boost::scoped_ptr<vw::SrcImageResource> src(vw::DiskImageResource::open(filename));
  return src->format();
}
