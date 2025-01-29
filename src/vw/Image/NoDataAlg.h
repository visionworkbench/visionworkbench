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

/// \file NoDataAlg.h
///
#ifndef __VW_IMAGE_NODATA_ALG_H__
#define __VW_IMAGE_NODATA_ALG_H__

#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PerPixelAccessorViews.h>

namespace vw {
  
// Function for highlighting spots of data
template<class PixelT>
class NotNoDataFunctor {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_nodata;
  typedef ChannelRange<channel_type> range_type;
public:
  NotNoDataFunctor( channel_type nodata ) : m_nodata(nodata) {}

  template <class Args> struct result {
    typedef channel_type type;
  };

  inline channel_type operator()( channel_type const& val ) const {
    return (val != m_nodata) && !std::isnan(val) ? range_type::max() : range_type::min();
  }
};

template <class ImageT, class NoDataT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type>  >
inline notnodata( ImageViewBase<ImageT> const& image, NoDataT nodata ) {
  typedef UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  func_type func( nodata );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}
  
}

#endif // __VW_IMAGE_NODATA_ALG_H__

