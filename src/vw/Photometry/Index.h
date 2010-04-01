// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Index.h

#ifndef __VW_PHOTOMETRY_INDEX_H__
#define __VW_PHOTOMETRY_INDEX_H__

#include <vector>
#include <string>

namespace vw {
namespace photometry {

  void index_images(std::vector<std::string> index_files,
                    std::vector<std::string> input_files);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_INDEX_H__
