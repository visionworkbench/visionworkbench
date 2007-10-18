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

#define DEFAULT_INTEREST_THRESHOLD (0.03)

namespace vw {
namespace ip {

  /// CRTP base class for thresholding interest points. Implementations
  /// must define operator(), which determines whether to keep or reject
  /// an interest point.
  template <class ImplT>
  struct InterestThresholdBase {
    /// Returns the derived implementation type.
    ImplT& impl() { return *static_cast<ImplT*>(this); }

    /// Returns the derived implementation type.
    ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

    // To be defined by the implementation.
    // template <class DataT>
    // inline bool operator() (InterestPoint const& pt, DataT const& data) const;
  };

  /// Simple thresholding class which keeps interest points above some
  /// minimum interest.
  template <int PeakT> struct InterestThreshold {};

  /// Specialization for thresholding maxima.
  template <> class InterestThreshold <IP_MAX> : public InterestThresholdBase<InterestThreshold<IP_MAX> > {
  private:
    float m_threshold;

  public:
    InterestThreshold(float min = DEFAULT_INTEREST_THRESHOLD) : m_threshold(min) {}

    template <class DataT>
    inline bool operator() (InterestPoint const& pt, DataT const& data) const {
      return (pt.interest > m_threshold);
    }
  };

  /// Specialization for thresholding minima.
  template <> class InterestThreshold <IP_MIN> : public InterestThresholdBase<InterestThreshold<IP_MIN> > {
  private:
    float m_threshold;

  public:
    InterestThreshold(float max = -DEFAULT_INTEREST_THRESHOLD) : m_threshold(max) {}

    template <class DataT>
    inline bool operator() (InterestPoint const& pt, DataT const& data) const {
      return (pt.interest < m_threshold);
    }
  };

  /// Specialization for thresholding both maxima and minima.
  template <> class InterestThreshold <IP_MINMAX> : public InterestThresholdBase<InterestThreshold<IP_MINMAX> > {
  private:
    float m_threshold;

  public:
    InterestThreshold(float threshold = DEFAULT_INTEREST_THRESHOLD) : m_threshold(threshold) {}

    template <class DataT>
    inline bool operator() (InterestPoint const& pt, DataT const& data) const {
      return (fabs(pt.interest) > m_threshold);
    }
  };

  // TODO: Better default thresholding options?  This is tricky to
  // decouple completely, as thresholding techniques can require data
  // not currently encapsulated by ImageInterestData, like curvatures
  // (from localization) and the Hessian matrix.

}} // namespace vw::ip

#endif //__INTERESTPOINT_THRESHOLD_H__
