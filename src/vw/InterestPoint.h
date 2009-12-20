// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file InterestPoint.h
///
/// A convenience header that includes the header files in vw/InterestPoint.
///
#ifndef __VW_INTERESTPOINT_H__
#define __VW_INTERESTPOINT_H__

// Data Types
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/InterestTraits.h>
#include <vw/InterestPoint/ImageOctave.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>
#include <vw/InterestPoint/IntegralImage.h>

// Operators
#include <vw/InterestPoint/InterestOperator.h>
#include <vw/InterestPoint/BoxFilter.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>

// Detection
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralDetector.h>

// Description
#include <vw/InterestPoint/LearnPCA.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/IntegralDescriptor.h>

// Matching
#include <vw/InterestPoint/Matcher.h>

#endif // __VW_INTERESTPOINT_H__

