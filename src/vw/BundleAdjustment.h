// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BundleAdjustment.h
///
/// Optimization classes for carrying out bundle adjustment of cameras.

#ifndef __VW_BUNDLE_ADJUSTMENT_H__
#define __VW_BUNDLE_ADJUSTMENT_H__

// Models
#include <vw/BundleAdjustment/ModelBase.h>

// Bundle Adjustment Implementations
#include <vw/BundleAdjustment/AdjustBase.h>
#include <vw/BundleAdjustment/AdjustRef.h>
#include <vw/BundleAdjustment/AdjustSparse.h>
#include <vw/BundleAdjustment/AdjustRobustRef.h>
#include <vw/BundleAdjustment/AdjustRobustSparse.h>

// Reporter
#include <vw/BundleAdjustment/BundleAdjustReport.h>

// Loading Utilities
#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/BundleAdjustment/ControlNetworkLoader.h>

#endif // __VW_BUNDLE_ADJUSTMENT_H__

