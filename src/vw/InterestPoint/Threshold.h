// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file Threshold.h
/// 
/// Functions for interest point thresholding.
/// 
#ifndef __INTERESTPOINT_THRESHOLD_H__
#define __INTERESTPOINT_THRESHOLD_H__

#include <vw/InterestPoint/InterestData.h>

namespace vw {
namespace ip {

  template <class ImplT>
  struct KeypointThresholdBase {
    /// Returns the derived implementation type.
    ImplT& impl() { return *static_cast<ImplT*>(this); }

    /// Returns the derived implementation type.
    ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

    template <class T>
    bool keep(const InterestPoint& pt,
              const ImageInterestData<T>& data) {
      return impl().keep(pt, data);
    }
  };

  template <class T>
  class InterestThreshold : public KeypointThresholdBase<InterestThreshold<T> > {
  public:
    T min_interest;

    InterestThreshold(T min = 0) : min_interest(min) {}

    bool keep(const InterestPoint& pt, const ImageInterestData<T>& data) {
      return (pt.interest > min_interest);
    }
  };

  // TODO: BrownLoweThreshold

}} // namespace vw::ip

#endif //__INTERESTPOINT_THRESHOLD_H__
